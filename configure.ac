AC_INIT(microsimulation, 1.7.6, [Mark Clements (firstname.lastname@ki.se)])

AC_SUBST([TPROCESS_IMPL], [0])
AC_CHECK_HEADER([setjmp.h], [AC_SUBST([TPROCESS_IMPL], [2])])
AC_CHECK_HEADER([ucontext.h], [AC_SUBST([TPROCESS_IMPL], [1])])

AC_CANONICAL_HOST

# turn off tprocess implementation for mac OS (deprecated)
case "${host_os}" in
     *darwin*)
	AC_SUBST([TPROCESS_IMPL], [0])
        ;;

     *)
	;;
esac

AC_OUTPUT([inst/include/tprocessconf.h])
