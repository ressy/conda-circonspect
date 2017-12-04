#!/usr/bin/env bash
set -Ee

envname=${1:-circonspect}

CIRCONSPECT_SRC="Circonspect-0.2.6"

function setup {
	# Setup Anaconda environment
	# Letting Circonspect try to handle the dependencies leads to a bunch
	# of compilation errors for me.  Just installing those via Anaconda
	# seems OK though.
	conda env list | grep ^"$envname " || \
		conda env create --file environment.yml --name "$envname" | tee conda_env_create.log
	source activate "$envname"
	# For some reason the Anaconda-provided perl doesn't automatically
	# figure out its library path.  We'll just force that automatically
	# during activation/deactivation of the environment.
	pushd "$CONDA_PREFIX"
	mkdir -p ./etc/conda/activate.d
	mkdir -p ./etc/conda/deactivate.d
	echo 'export PERL5LIB_ORIG="$PERL5LIB"; export PERL5LIB="$CONDA_PREFIX/lib/perl5"' > ./etc/conda/activate.d/perl.sh
	echo 'export PERL5LIB="$PERL5LIB_ORIG"; unset PERL5LIB_ORIG' > ./etc/conda/deactivate.d/perl.sh
	popd

	# Setup Circonspect
	pushd "$CIRCONSPECT_SRC"
	[[ -e Makefile ]] && make clean
	# We'll set the install prefix to the Anaconda environment root
	yes | perl Makefile.PL INSTALL_BASE="$CONDA_PREFIX" | tee ../perl_Makefile.PL.log
	make | tee ../make.log
	make test | tee ../make_test.log
	make install | tee ../make_install.log
	popd
}

function catch {
	echo ""
	echo "Error during setup"
	exit 1
}

function post_setup {
echo
echo "To use:"
echo "$ source activate $envname"
echo "$ Circonspect -help"
}

trap catch ERR
setup
post_setup
