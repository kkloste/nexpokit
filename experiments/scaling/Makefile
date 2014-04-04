
CXXFLAGS += -std=c++0x -O3 -Wall
CC := $(CXX)

CXX=g++-4.4

all :  forest_fire_model

forest_fire_model: forest_fire_model.o sf_util.o


.PHONE : clean

clean :
	rm -rf sf_util.o \
      copying_model copying_model.o \
      preferential_attachment_model.o preferential_attachment_model \
      forest_fire_model forest_fire_model.o


