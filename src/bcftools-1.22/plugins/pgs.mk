# The pgs plugin requires CHOLMOD library for sparse matrix operations
# Only build if CHOLMOD is available
HAVE_CHOLMOD=yes
CHOLMOD_LIBS=-L/home/sounkoutoure/Projects/RBCFLib/src/SuiteSparse/install/lib -lcholmod -lccolamd -lcolamd -lcamd -lamd -lsuitesparseconfig
CHOLMOD_CPPFLAGS=-I/home/sounkoutoure/Projects/RBCFLib/src/SuiteSparse/install/include -I/home/sounkoutoure/Projects/RBCFLib/src/SuiteSparse/install/include/suitesparse
ifeq ($(HAVE_CHOLMOD),yes)
    plugins/pgs.so: plugins/pgs.c
	    $(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(ALL_CPPFLAGS) $(CHOLMOD_CPPFLAGS) $(EXTRA_CPPFLAGS) $(LDFLAGS) -o $@ version.c $< $(PLUGIN_LIBS) $(LIBS) $(CHOLMOD_LIBS)
    plugins/pgs.dll: plugins/pgs.c
	    $(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(ALL_CPPFLAGS) $(CHOLMOD_CPPFLAGS) $(EXTRA_CPPFLAGS) $(LDFLAGS) -o $@ version.c $< $(PLUGIN_LIBS) $(LIBS) $(CHOLMOD_LIBS)
else
    # Create a dummy target if CHOLMOD is not available
    plugins/pgs.so:
	    @echo "Warning: pgs plugin requires CHOLMOD library, skipping build"
	    @touch $@
    plugins/pgs.dll:
	    @echo "Warning: pgs plugin requires CHOLMOD library, skipping build"
	    @touch $@
		
endif