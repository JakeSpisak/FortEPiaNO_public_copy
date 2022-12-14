default: fortepiano

all: fortepiano tests pythonwrapper

fortepiano: BUILD_DIR ?= build
readnodes: BUILD_DIR ?= build
pythonwrapper: BUILD_DIR ?= build
preparenodes: BUILD_DIR ?= build
prepareinterp: BUILD_DIR ?= build
tests: BUILD_DIR ?= buildtest

fortepiano: directories
	cd ./sources && make fortepiano BUILD_DIR=../$(BUILD_DIR)

readnodes: directories
	cd ./sources && make readnodes BUILD_DIR=../$(BUILD_DIR)

pythonwrapper: directories
	cd ./sources && make pythonwrapper BUILD_DIR=../$(BUILD_DIR)

preparenodes: directories
	cd ./sources && make preparenodes BUILD_DIR=../$(BUILD_DIR)

prepareinterp: directories
	cd ./sources && make prepareinterp BUILD_DIR=../$(BUILD_DIR)

tests: directories
	cd ./sources && make tests BUILD_DIR=../$(BUILD_DIR)

directories:
	mkdir -p bin/ log/ $(BUILD_DIR) GL_nodes/ interpolations/

clean: 
	rm -rf build*/

cleanwrapper: clean
	rm -rf python/fortepianoWrapper*.so sources/fortepianoWrapper*.so

cleanall: clean cleanwrapper
	rm -rf bin/*

installpython:
	python setup.py install --user
