#!/bin/bash -x

# exit when any command fails
set -e

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <version> <install_dir> <install_dep_dir>"
    exit 1
fi

VERSION=$1
INSTALL_DIR=$2
INSTALL_DEP_DIR=$3
ROOT=`pwd`

if [[ "$OSTYPE" == "darwin"* ]]; then
    PLATFORM="mac"
else
    PLATFORM="linux"
fi

echo "VERSION = $VERSION"
echo "PR_NUMBER = $PR_NUMBER"
echo "GITHUB_REF = $GITHUB_REF"
echo "PLATFORM = $PLATFORM"

if [[ "$PR_NUMBER" != "" ]]; then
    VERSION="ShapeWorks-PR-${PR_NUMBER}-${PLATFORM}"
else
    if [[ "$VERSION" == "tag" ]]; then
	VERSION="ShapeWorks-$(git describe --tags)-${PLATFORM}"
    fi
fi

# Special case for when we are on the master branch (dev releases)
BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [[ "$BRANCH" == "master" ]]; then
    NUMBER=$(git rev-list start_dev_releases..HEAD --count)
    HASH=$(git rev-parse --short HEAD)
    VERSION="ShapeWorks-dev-${NUMBER}-${HASH}-${PLATFORM}"
fi

echo "Version: $VERSION"

rm -rf "package/$VERSION"

mkdir -p "package/$VERSION"

# Build python package tarballs
# Pip can't install these otherwise from a read-only area like /Applications
for package in DataAugmentationUtilsPackage DatasetUtilsPackage DeepSSMUtilsPackage DocumentationUtilsPackage ShapeCohortGenPackage shapeworks ; do
    cd Python
    tar czvf ${package}.tar.gz $package
    cd ..
done

BASE_LIB=${INSTALL_DEP_DIR}/lib
cp -a $INSTALL_DIR/* "package/${VERSION}"
cp -a Examples "package/${VERSION}"
cp -a Python "package/${VERSION}"
cp -a Installation "package/${VERSION}"
cp install_shapeworks.sh package/${VERSION}
cp python_requirements.txt package/${VERSION}
cp docs/about/release-notes.md package/${VERSION}

if [[ "$OSTYPE" == "darwin"* ]]; then
    cp -a $INSTALL_DEP_DIR/lib/*.dylib "package/${VERSION}/bin/ShapeWorksStudio.app/Contents/Frameworks"
    cp docs/users/Mac_README.txt package/${VERSION}/README.txt
    if [ $? -ne 0 ]; then
	echo "Failed to copy Mac package README"
	exit 1
    fi
else
    cp -a $INSTALL_DEP_DIR/lib "package/${VERSION}"
    cp docs/users/Linux_README.txt package/${VERSION}/README.txt
    if [ $? -ne 0 ]; then
	echo "Failed to copy Linux package README"
	exit 1
    fi
fi

cd "package/${VERSION}"
rm -rf include share v3p plugins libigl

if [[ "$OSTYPE" == "darwin"* ]]; then
    # Mac OSX
    cd bin
    install_name_tool -add_rpath @executable_path/../Frameworks ShapeWorksStudio.app/Contents/MacOS/ShapeWorksStudio || echo ok
    QT_LIB_LOCATION="@executable_path/ShapeWorksStudio.app/Contents/Frameworks"
    QT_LOADER_LIB_LOCATION="@loader_path/ShapeWorksStudio.app/Contents/Frameworks"

    # copy platform plugins for Studio
    cp -a ShapeWorksStudio.app/Contents/PlugIns .

    for i in *.so ; do
	install_name_tool -add_rpath "@loader_path/../lib" $i || echo ok
	install_name_tool -add_rpath $QT_LOADER_LIB_LOCATION $i || echo ok
    done

    for i in * ; do
	install_name_tool -add_rpath $QT_LIB_LOCATION $i || echo ok
    done

    # Copy libraries from anaconda
    conda_libs="libpython libboost_filesystem"
    for clib in $conda_libs; do
        cp ${CONDA_PREFIX}/lib/${clib}* ShapeWorksStudio.app/Contents/Frameworks
    done


    # # Fix transitive loaded libs
    # for i in ShapeWorksStudio.app/Contents/Frameworks/*.dylib ; do
    # 	install_name_tool -change ${BASE_LIB}/libitkgdcmopenjp2-5.2.1.dylib @rpath/libitkgdcmopenjp2-5.2.1.dylib $i
    # done
    # install_name_tool -id @rpath/libitkgdcmopenjp2-5.2.1.dylib ShapeWorksStudio.app/Contents/Frameworks/libitkgdcmopenjp2-5.2.1.dylib

    cd ..

    # remove static libs
    pwd
    rm lib/*.a

else
    # Copy libraries from anaconda
    conda_libs="libboost_iostreams libboost_filesystem libbz2 liblzma libtbb libHalf libpython libz libspd"
    for clib in $conda_libs; do
        cp ${CONDA_PREFIX}/lib/${clib}* lib
    done

    # remove static libs
    rm lib/*.a
    
    cd bin
    linuxdeployqt ShapeWorksStudio -verbose=2
    cd ..
    
    rm lib/libxcb* lib/libX* lib/libfont* lib/libfreetype*
    rm -rf geometry-central doc
fi

# Run auto-documentation
cd $ROOT
PATH=$ROOT/package/${VERSION}/bin:$PATH
# check that 'shapeworks -h' is working
shapeworks -h
if [ $? -eq 0 ]; then
    echo "shapeworks -h is working"
else
    echo "shapeworks -h is not working"
    exit 1
fi
python Python/RunShapeWorksAutoDoc.py --md_filename docs/tools/ShapeWorksCommands.md
mkdocs build
mv site Documentation
cp -a Documentation "${ROOT}/package/${VERSION}"

mkdir ${ROOT}/artifacts
cd ${ROOT}/package
if [[ "$OSTYPE" == "darwin"* ]]; then
    zip -y -r ${ROOT}/artifacts/${VERSION}.zip ${VERSION}
    if [ $? -ne 0 ]; then
	echo "Failed to zip artifact"
	exit 1
    fi
else
    tar czvf ${ROOT}/artifacts/${VERSION}.tar.gz ${VERSION}
    if [ $? -ne 0 ]; then
	echo "Failed to tar artifact"
	exit 1
    fi
fi

# Additionally on Mac, create an installer
if [[ "$OSTYPE" == "darwin"* ]]; then
    cp ${ROOT}/docs/users/Mac_README.txt ${VERSION}/README.txt
    pkgbuild --quiet --analyze --root ${VERSION} ShapeWorks.plist
    plutil -replace BundleIsRelocatable -bool NO ShapeWorks.plist
    pkgbuild --component-plist ShapeWorks.plist --install-location /Applications/ShapeWorks --root ${VERSION} --identifier edu.utah.sci.shapeworks ${ROOT}/artifacts/${VERSION}.pkg --scripts ${ROOT}/Support/osxscripts
fi

cd $ROOT


