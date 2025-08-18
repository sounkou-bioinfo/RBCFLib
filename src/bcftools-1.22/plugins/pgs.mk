# The pgs plugin requires CHOLMOD library for sparse matrix operations
# Only build if CHOLMOD is available
ifdef HAVE_CHOLMOD
    ifndef CHOLMOD_LIBS
        CHOLMOD_LIBS = -lcholmod
    endif
    plugins/pgs.so: plugins/pgs.c
	    $(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(ALL_CPPFLAGS) $(EXTRA_CPPFLAGS) $(LDFLAGS) -o $@ version.c $< $(PLUGIN_LIBS) $(LIBS) $(CHOLMOD_LIBS)
else
    # Create a dummy target if CHOLMOD is not available
    plugins/pgs.so:
	    @echo "Warning: pgs plugin requires CHOLMOD library, skipping build"
	    @touch $@
endif