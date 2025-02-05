CXX=g++
CXX_COMMON_LIB=-Wall -Werror -Wpedantic -pipe -O2 -fPIC -c
CXX_COMMON_EXE=-Wall -Wpedantic -pipe -O2 -g -ldl

JSON_INCLUDE=-I /usr/include/jsoncpp
JSON_LIB=-ljsoncpp

CPP_TOOLS_LIB_DIR=./CppTools/lib
CPP_TOOLS_SRC_DIR=./CppTools/src
CPP_TOOLS_INCLUDE=-I ./CppTools/include
CPP_TOOLS_LIB=-L $(CPP_TOOLS_LIB_DIR) -Wl,-rpath,$(CPP_TOOLS_LIB_DIR) -lBox -lIOTools -lStrTools -lErrorHandler 

ROOT_TOOLS_LIB_DIR=./ROOTTools/lib
ROOT_TOOLS_SRC_DIR=./ROOTTools/src
ROOT_TOOLS_INCLUDE=-I ./ROOTTools/include
ROOT_TOOLS_LIB=-L $(ROOT_TOOLS_LIB_DIR) -Wl,-rpath,$(ROOT_TOOLS_LIB_DIR) -lTCanvasTools -lTF1Tools

PBAR_LIB_DIR=./ProgressBar/lib
PBAR_SRC_DIR=./ProgressBar/src
PBAR_INCLUDE=-I ./ProgressBar/include
PBAR_LIB=-L $(PBAR_LIB_DIR) -Wl,-rpath,$(PBAR_LIB_DIR) -lPBar

ROOT_INCLUDE=-I${ROOT_PATH}/include
ROOT_CONFIG=${ROOT_PATH}/bin/root-config

CAL_PHENIX_INCLUDE=-I ./include
CAL_PHENIX_LIB=-L ./lib -Wl,-rpath,./lib -lInputReader

ALL_INCLUDE=$(JSON_INCLUDE) $(CPP_TOOLS_INCLUDE) $(ROOT_INCLUDE) $(PBAR_INCLUDE) $(ROOT_TOOLS_INCLUDE) $(CAL_PHENIX_INCLUDE)
ALL_LIB=$(JSON_LIB) $(CPP_TOOLS_LIB) `$(ROOT_CONFIG) --glibs` $(PBAR_LIB) $(ROOT_TOOLS_LIB) $(CAL_PHENIX_LIB)
