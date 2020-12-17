ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = foreign
SUBDIRS = src tests
AM_CXXFLAGS = $(MASA_CXXFLAGS) $(GRVY_CFLAGS) $(HDF5_CFLAGS) 
LIBS = $(MASA_LIBS) $(GRVY_LIBS) $(HDF5_LIBS)
MKDIR_P = mkdir -p
if CODE_COVERAGE_ENABLED

lcov_dir = $(top_builddir)/coverage/lcov

lcov-report:
	$(MKDIR_P) $(lcov_dir)
	$(top_srcdir)/src/lcov/bin/lcov --compat-libtool --directory . --capture --output-file $(lcov_dir)/lcov.info
	$(top_srcdir)/src/lcov/bin/lcov --list-full-path -l $(lcov_dir)/lcov.info | grep -v "`cd -P $(top_srcdir)/src && pwd`" | cut -d\| -f1 > $(lcov_dir)/remove
	$(top_srcdir)/src/lcov/bin/lcov -q -r $(lcov_dir)/lcov.info `cat $(lcov_dir)/remove` > $(lcov_dir)/lcov.cleaned.info
	@rm $(lcov_dir)/remove
	@mv $(lcov_dir)/lcov.cleaned.info $(lcov_dir)/lcov.info
	$(top_srcdir)/src/lcov/bin/genhtml -t "Project 1" -o $(lcov_dir) $(lcov_dir)/lcov.info

lcov-reset:
	@rm -rf $(lcov_dir)
	@find . -name "*.gcda" -exec rm {} \;
	$(top_srcdir)/src/lcov/bin/lcov --directory . --zerocounters

coverage: tests lcov-reset check lcov-report

endif
