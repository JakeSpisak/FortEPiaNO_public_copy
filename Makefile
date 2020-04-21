default: fortepiano

all: fortepiano tests

fortepiano: BUILD_DIR ?= build
tests: BUILD_DIR ?= buildtest

fortepiano: directories
	cd ./sources && make fortepiano BUILD_DIR=../$(BUILD_DIR)

tests: directories
	cd ./sources && make tests BUILD_DIR=../$(BUILD_DIR)

directories:
	mkdir -p bin/ log/ $(BUILD_DIR)

clean: 
	rm -rf bin/* build*/

installpython:
	python setup.py install --user
