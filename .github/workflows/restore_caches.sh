#!/bin/bash -x

echo "#############################"
echo "# Restore Caches            #"
echo "#############################"

. $GITHUB_WORKSPACE/.github/workflows/common.sh

# based on root folder
cd /

if [[ "$USE_CCACHE" == "ON" ]]; then
    # Restore ccache
    scp runner@${CACHE_HOST}:github/${CCACHE_FILE} /tmp
    if [ -f ${CCACHE_FILE} ] ; then
	echo "ccache file was found"
	tar --use-compress-program=pigz -xf /tmp/${CCACHE_FILE}
	rm /tmp/$CCACHE_FILE
    fi
fi

# Restore conda installs
scp runner@${CACHE_HOST}:github/${CONDA_FILE} /tmp
if [ -f ${CONDA_FILE} ] ; then
    echo "Conda Cache file was found"
    tar --use-compress-program=pigz -xf /tmp/$CONDA_FILE
    rm /tmp/$CONDA_FILE
fi

# Restore dependencies
scp runner@${CACHE_HOST}:github/${DEP_FILE} /tmp
if [ -f ${DEP_FILE} ] ; then
    echo "Dependency Cache file was found"
    tar --use-compress-program=pigz -xf /tmp/$DEP_FILE
    rm /tmp/$DEP_FILE
fi

which ccache
ccache -p
ccache -s
