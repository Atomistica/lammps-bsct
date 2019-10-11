# Install/unInstall package files in LAMMPS

if (test $1 = 2) then

    for i in anderson_mixer.h anderson_mixer.cpp fix_bsct.h fix_bsct.cpp fix_bsct_parameters_t.h pair_lj_charmmfsw_coul_charmmfsh_bsct.h pair_lj_charmmfsw_coul_charmmfsh_bsct.cpp pair_lj_charmmfsw_coul_long_bsct.h pair_lj_charmmfsw_coul_long_bsct.cpp pppm_bsct.h pppm_bsct.cpp ; do
        diff -q ../$i $i
    done

elif (test $1 = 1) then

    for i in anderson_mixer.h anderson_mixer.cpp fix_bsct.h fix_bsct.cpp fix_bsct_parameters_t.h pair_lj_charmmfsw_coul_charmmfsh_bsct.h pair_lj_charmmfsw_coul_charmmfsh_bsct.cpp pair_lj_charmmfsw_coul_long_bsct.h pair_lj_charmmfsw_coul_long_bsct.cpp pppm_bsct.h pppm_bsct.cpp ; do
	cp -p $i ..
    done

    if (test -e ../Makefile.package) then
	sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(user-bsct_SYSINC) |' ../Makefile.package
	sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(user-bsct_SYSLIB) |' ../Makefile.package
	sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(user-bsct_SYSPATH) |' ../Makefile.package
    fi

elif (test $1 = 0) then

    for i in anderson_mixer.h anderson_mixer.cpp fix_bsct.h fix_bsct.cpp fix_bsct_parameters_t.h pair_lj_charmmfsw_coul_charmmsh_bsct.h pair_lj_charmmfsw_coul_charmmfsh_bsct.cpp pair_lj_charmmfsw_coul_long_bsct.h pair_lj_charmmfsw_coul_long_bsct.cpp pppm_bsct.h pppm_bsct.cpp ; do
	rm -f ../$i
    done

    if (test -e ../Makefile.package) then
	sed -i -e 's/[^ \t]*bsct[^ \t]* //' ../Makefile.package
    fi

fi
