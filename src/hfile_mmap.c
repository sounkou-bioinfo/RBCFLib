/*  hfile_mmap.c -- Memory-mapped local file backend for low-level file streams.

    Copyright (C) 2016 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifdef _WIN32
#include <windows.h>
#include <io.h>
#include <fcntl.h>
#include <sys/stat.h>
#include "./bcftools-1.22/htslib-1.22/hfile_internal.h"
#endif

#ifndef _WIN32
#include <errno.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "./bcftools-1.22/htslib-1.22/hfile_internal.h"

typedef struct {
    hFILE base;
    char *buffer;
    size_t length, pos;
    int fd;
} hFILE_mmap;

static ssize_t mmap_read(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_mmap *fp = (hFILE_mmap *) fpv;
    size_t avail = fp->length - fp->pos;
    if (nbytes > avail) nbytes = avail;
    memcpy(buffer, fp->buffer + fp->pos, nbytes);
    fp->pos += nbytes;
    return nbytes;
}

static ssize_t mmap_write(hFILE *fpv, const void *buffer, size_t nbytes)
{
    hFILE_mmap *fp = (hFILE_mmap *) fpv;
    size_t avail = fp->length - fp->pos;
    if (nbytes > avail) nbytes = avail;
    memcpy(fp->buffer + fp->pos, buffer, nbytes);
    fp->pos += nbytes;
    return nbytes;
}

static off_t mmap_seek(hFILE *fpv, off_t offset, int whence)
{
    hFILE_mmap *fp = (hFILE_mmap *) fpv;
    size_t absoffset = (offset >= 0)? offset : -offset;
    size_t origin;
    //fprintf(stderr, "[mmap] seek %ld bytes %s\n", absoffset, (offset < 0)? "backward" : "forward");

    switch (whence) {
    case SEEK_SET: origin = 0; break;
    case SEEK_CUR: origin = fp->pos; break;
    case SEEK_END: origin = fp->length; break;
    default: errno = EINVAL; return -1;
    }

    if ((offset  < 0 && absoffset > origin) ||
        (offset >= 0 && absoffset > fp->length - origin)) {
        errno = EINVAL;
        return -1;
    }

    fp->pos = origin + offset;
    return fp->pos;
}

static int mmap_close(hFILE *fpv)
{
    hFILE_mmap *fp = (hFILE_mmap *) fpv;
    int ret = 0;
    if (munmap(fp->buffer, fp->length) < 0) ret = -1;
    if (close(fp->fd) < 0) ret = -1;
    return ret;
}

static const struct hFILE_backend mmap_backend =
{
    mmap_read, mmap_write, mmap_seek, NULL, mmap_close
};

static hFILE *hopen_mmap(const char *filename, const char *modestr)
{
    int mode = hfile_oflags(modestr);
    struct stat st;
    int fd = -1;
    void *data = MAP_FAILED;
    hFILE_mmap *fp = NULL;
    int prot, save;

    if (strncmp(filename, "mmap://localhost/", 17) == 0) filename += 16;
    else if (strncmp(filename, "mmap:///", 8) == 0) filename += 7;
    else if (strncmp(filename, "mmap:", 5) == 0) filename += 5;

    fd = open(filename, mode, 0666);
    if (fd < 0) goto error;
    if (fstat(fd, &st) < 0) goto error;

    switch (mode & O_ACCMODE) {
    case O_RDONLY: prot = PROT_READ;  break;
    case O_WRONLY: prot = PROT_WRITE; break;
    case O_RDWR:   prot = PROT_READ|PROT_WRITE; break;
    default:       prot = PROT_NONE;  break;
    }

    data = mmap(NULL, st.st_size, prot, MAP_SHARED, fd, 0);
//    fprintf(stderr, "[mmap] mapped %s at %p, length=%zu bytes\n",
  //      filename, data, (size_t) st.st_size);
    if (data == MAP_FAILED) goto error;

    fp = (hFILE_mmap *) hfile_init(sizeof (hFILE_mmap), modestr, st.st_blksize);
    if (fp == NULL) goto error;

    fp->fd = fd;
    fp->buffer = data;
    fp->length = st.st_size;
    fp->pos = 0;
    fp->base.backend = &mmap_backend;
    return &fp->base;

error:
    save = errno;
    if (fp) hfile_destroy((hFILE *) fp);
    if (data != MAP_FAILED) (void) munmap(data, st.st_size);
    if (fd >= 0) (void) close(fd);
    errno = save;
    return NULL;
}

int hfile_plugin_init(struct hFILE_plugin *self)
{
    static const struct hFILE_scheme_handler handler =
        { hopen_mmap, hfile_always_local, "mmap", 10 };

    if (self) self->name = "mmap";
    hfile_add_scheme_handler("mmap", &handler);
    return 0;
}

// Wrapper function for easier calling
int hfile_plugin_init_mmap(void) {
    return hfile_plugin_init(NULL);
}
#endif

#ifdef _WIN32
// Windows implementation using CreateFileMapping and MapViewOfFile

typedef struct {
    hFILE base;
    char *buffer;
    size_t length, pos;
    HANDLE file_handle;
    HANDLE mapping_handle;
} hFILE_mmap_win;

static ssize_t mmap_read_win(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_mmap_win *fp = (hFILE_mmap_win *) fpv;
    size_t avail = fp->length - fp->pos;
    if (nbytes > avail) nbytes = avail;
    memcpy(buffer, fp->buffer + fp->pos, nbytes);
    fp->pos += nbytes;
    return nbytes;
}

static ssize_t mmap_write_win(hFILE *fpv, const void *buffer, size_t nbytes)
{
    hFILE_mmap_win *fp = (hFILE_mmap_win *) fpv;
    size_t avail = fp->length - fp->pos;
    if (nbytes > avail) nbytes = avail;
    memcpy(fp->buffer + fp->pos, buffer, nbytes);
    fp->pos += nbytes;
    return nbytes;
}

static off_t mmap_seek_win(hFILE *fpv, off_t offset, int whence)
{
    hFILE_mmap_win *fp = (hFILE_mmap_win *) fpv;
    size_t absoffset = (offset >= 0)? offset : -offset;
    size_t origin;

    switch (whence) {
    case SEEK_SET: origin = 0; break;
    case SEEK_CUR: origin = fp->pos; break;
    case SEEK_END: origin = fp->length; break;
    default: 
        SetLastError(ERROR_INVALID_PARAMETER);
        return -1;
    }

    if ((offset  < 0 && absoffset > origin) ||
        (offset >= 0 && absoffset > fp->length - origin)) {
        SetLastError(ERROR_INVALID_PARAMETER);
        return -1;
    }

    fp->pos = origin + offset;
    return fp->pos;
}

static int mmap_close_win(hFILE *fpv)
{
    hFILE_mmap_win *fp = (hFILE_mmap_win *) fpv;
    int ret = 0;
    
    if (fp->buffer) {
        if (!UnmapViewOfFile(fp->buffer)) ret = -1;
    }
    if (fp->mapping_handle != NULL && fp->mapping_handle != INVALID_HANDLE_VALUE) {
        if (!CloseHandle(fp->mapping_handle)) ret = -1;
    }
    if (fp->file_handle != NULL && fp->file_handle != INVALID_HANDLE_VALUE) {
        if (!CloseHandle(fp->file_handle)) ret = -1;
    }
    return ret;
}

static const struct hFILE_backend mmap_backend_win =
{
    mmap_read_win, mmap_write_win, mmap_seek_win, NULL, mmap_close_win
};

static hFILE *hopen_mmap_win(const char *filename, const char *modestr)
{
    int mode = hfile_oflags(modestr);
    HANDLE file_handle = INVALID_HANDLE_VALUE;
    HANDLE mapping_handle = INVALID_HANDLE_VALUE;
    void *data = NULL;
    hFILE_mmap_win *fp = NULL;
    DWORD access, protect, map_access;
    LARGE_INTEGER file_size;

    // Strip mmap:// prefixes
    if (strncmp(filename, "mmap://localhost/", 17) == 0) filename += 16;
    else if (strncmp(filename, "mmap:///", 8) == 0) filename += 7;
    else if (strncmp(filename, "mmap:", 5) == 0) filename += 5;

    // Convert Unix-style file access modes to Windows
    switch (mode & O_ACCMODE) {
    case O_RDONLY: 
        access = GENERIC_READ;
        protect = PAGE_READONLY;
        map_access = FILE_MAP_READ;
        break;
    case O_WRONLY:
    case O_RDWR:   
        access = GENERIC_READ | GENERIC_WRITE;
        protect = PAGE_READWRITE;
        map_access = FILE_MAP_WRITE;
        break;
    default:       
        SetLastError(ERROR_INVALID_PARAMETER);
        goto error;
    }

    // Open the file
    file_handle = CreateFile(filename, access, FILE_SHARE_READ, NULL, 
                            OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (file_handle == INVALID_HANDLE_VALUE) goto error;

    // Get file size
    if (!GetFileSizeEx(file_handle, &file_size)) goto error;
    if (file_size.QuadPart == 0) {
        // Handle empty files
        fp = (hFILE_mmap_win *) hfile_init(sizeof(hFILE_mmap_win), modestr, 4096);
        if (fp == NULL) goto error;
        
        fp->file_handle = file_handle;
        fp->mapping_handle = NULL;
        fp->buffer = NULL;
        fp->length = 0;
        fp->pos = 0;
        fp->base.backend = &mmap_backend_win;
        return &fp->base;
    }

    // Create file mapping
    mapping_handle = CreateFileMapping(file_handle, NULL, protect, 
                                      file_size.HighPart, file_size.LowPart, NULL);
    if (mapping_handle == NULL) goto error;

    // Map view of file
    data = MapViewOfFile(mapping_handle, map_access, 0, 0, 0);
    if (data == NULL) goto error;

    fp = (hFILE_mmap_win *) hfile_init(sizeof(hFILE_mmap_win), modestr, 4096);
    if (fp == NULL) goto error;

    fp->file_handle = file_handle;
    fp->mapping_handle = mapping_handle;
    fp->buffer = (char*)data;
    fp->length = (size_t)file_size.QuadPart;
    fp->pos = 0;
    fp->base.backend = &mmap_backend_win;
    return &fp->base;

error:
    if (fp) hfile_destroy((hFILE *) fp);
    if (data) UnmapViewOfFile(data);
    if (mapping_handle != NULL && mapping_handle != INVALID_HANDLE_VALUE) 
        CloseHandle(mapping_handle);
    if (file_handle != NULL && file_handle != INVALID_HANDLE_VALUE) 
        CloseHandle(file_handle);
    return NULL;
}

int hfile_plugin_init(struct hFILE_plugin *self)
{
    static const struct hFILE_scheme_handler handler =
        { hopen_mmap_win, hfile_always_local, "mmap", 10 };

    if (self) self->name = "mmap";
    hfile_add_scheme_handler("mmap", &handler);
    return 0;
}

// Wrapper function for easier calling
int hfile_plugin_init_mmap(void) {
    return hfile_plugin_init(NULL);
}
#endif