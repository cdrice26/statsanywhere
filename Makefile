# Makefile for statsanywhere library (build from src/ into dist/)

CC ?= cc
AR ?= ar
RANLIB ?= ranlib
CFLAGS ?= -O2 -Wall -Wextra -std=c11
PICFLAGS := -fPIC
LDFLAGS ?=

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	SHLIB_EXT := .dylib
	SHLIB_LDFLAGS := -dynamiclib -install_name @rpath/libstatsanywhere$(SHLIB_EXT)
else
	SHLIB_EXT := .so
	SHLIB_LDFLAGS := -shared
endif

SRC_DIR := src
DIST_DIR := dist
OBJ_DIR := $(DIST_DIR)/obj

# All .c sources in src/ except main.c (main.c is for local testing/unity and excluded)
SRCS := $(filter-out $(SRC_DIR)/main.c,$(wildcard $(SRC_DIR)/*.c))
OBJS := $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRCS))

LIBNAME := statsanywhere
STATIC_LIB := $(DIST_DIR)/lib$(LIBNAME).a
SHARED_LIB := $(DIST_DIR)/lib$(LIBNAME)$(SHLIB_EXT)

.PHONY: all lib static shared clean distclean info

all: lib

lib: static shared

static: $(STATIC_LIB)

$(STATIC_LIB): CFLAGS += $(PICFLAGS) -I$(SRC_DIR)
$(STATIC_LIB): $(OBJS)
	@mkdir -p $(DIST_DIR)
	$(AR) rcs $@ $(OBJS)
	@$(RANLIB) $@ 2>/dev/null || true

shared: $(SHARED_LIB)

$(SHARED_LIB): CFLAGS += $(PICFLAGS) -I$(SRC_DIR)
$(SHARED_LIB): $(OBJS)
	@mkdir -p $(DIST_DIR)
	$(CC) $(SHLIB_LDFLAGS) -o $@ $(OBJS) $(LDFLAGS)

# Compile objects into $(OBJ_DIR)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -I$(SRC_DIR) -c $< -o $@

# Convenience rules
clean:
	rm -rf $(OBJ_DIR)
	rm -f $(DIST_DIR)/lib$(LIBNAME).a $(DIST_DIR)/lib$(LIBNAME)$(SHLIB_EXT)

distclean: clean
	rm -rf $(DIST_DIR)

info:
	@echo "Source dir: $(SRC_DIR)"
	@echo "Distribution dir: $(DIST_DIR)"
	@echo "Sources: $(SRCS)"
	@echo "Objects: $(OBJS)"
	@echo "Static library: $(STATIC_LIB)"
	@echo "Shared library: $(SHARED_LIB)"
	@echo "CC=$(CC) CFLAGS=$(CFLAGS)"
