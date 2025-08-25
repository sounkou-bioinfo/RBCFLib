HTSLIB_static_LDFLAGS = -Wl,-Bsymbolic-functions -flto=auto -ffat-lto-objects -Wl,-z,relro
HTSLIB_static_LIBS = -lpthread -lz -lm -lbz2 -llzma -ldeflate -lcurl -lcrypto
