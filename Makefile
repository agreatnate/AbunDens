# Makefile for Number Theory Density Calculation Programs
# Builds both abundant and covering number density calculators using shared framework
# Requires: GMP, MPFR, OpenMP

# ============================================================================
# Configuration
# ============================================================================

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -flto -DNDEBUG -fopenmp -frounding-math
DEBUGFLAGS = -std=c++17 -O0 -g -Wall -Wextra -fsanitize=address -fopenmp -frounding-math
INCLUDES = -I.
LIBS = -lgmp -lgmpxx -lmpfr

# Directories
SRCDIR = .
OBJDIR = obj
BINDIR = bin
DEPDIR = deps

# Create directories if they don't exist
$(shell mkdir -p $(OBJDIR) $(BINDIR) $(DEPDIR))

# ============================================================================
# Source Files
# ============================================================================

# Shared framework sources
SHARED_SOURCES = DensityMathFramework.cpp

# Abundant number sources
ABUND_SOURCES = AbunDens.cpp

# Covering number sources
COV_SOURCES = CovDens.cpp

# All sources
ALL_SOURCES = $(SHARED_SOURCES) $(ABUND_SOURCES) $(COV_SOURCES)

# ============================================================================
# Object Files
# ============================================================================

# Shared framework objects
SHARED_OBJECTS = $(SHARED_SOURCES:%.cpp=$(OBJDIR)/%.o)

# Abundant number objects
ABUND_OBJECTS = $(ABUND_SOURCES:%.cpp=$(OBJDIR)/%.o)

# Covering number objects
COV_OBJECTS = $(COV_SOURCES:%.cpp=$(OBJDIR)/%.o)

# All objects
ALL_OBJECTS = $(ALL_SOURCES:%.cpp=$(OBJDIR)/%.o)

# ============================================================================
# Dependency Files
# ============================================================================

DEPS = $(ALL_SOURCES:%.cpp=$(DEPDIR)/%.d)

# ============================================================================
# Targets
# ============================================================================

# Default target - build both programs
.PHONY: all
all: $(BINDIR)/abundens $(BINDIR)/covdens

# Individual program targets
.PHONY: abundant covering
abundant: $(BINDIR)/abundens
covering: $(BINDIR)/covdens

# Debug versions
.PHONY: debug debug-abundant debug-covering
debug: debug-abundant debug-covering
debug-abundant: $(BINDIR)/abundens-debug
debug-covering: $(BINDIR)/covdens-debug

# ============================================================================
# Program Linking Rules
# ============================================================================

# Abundant number density calculator
$(BINDIR)/abundens: $(SHARED_OBJECTS) $(ABUND_OBJECTS)
	@echo "Linking abundant number density calculator..."
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)
	@echo "Built: $@"

# Covering number density calculator
$(BINDIR)/covdens: $(SHARED_OBJECTS) $(COV_OBJECTS)
	@echo "Linking covering number density calculator..."
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)
	@echo "Built: $@"

# Debug versions
$(BINDIR)/abundens-debug: $(SHARED_OBJECTS:.o=-debug.o) $(ABUND_OBJECTS:.o=-debug.o)
	@echo "Linking debug abundant number calculator..."
	$(CXX) $(DEBUGFLAGS) -o $@ $^ $(LIBS)
	@echo "Built: $@"

$(BINDIR)/covdens-debug: $(SHARED_OBJECTS:.o=-debug.o) $(COV_OBJECTS:.o=-debug.o)
	@echo "Linking debug covering number calculator..."
	$(CXX) $(DEBUGFLAGS) -o $@ $^ $(LIBS)
	@echo "Built: $@"

# ============================================================================
# Compilation Rules
# ============================================================================

# Standard compilation rule with dependency generation
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -MF $(DEPDIR)/$*.d -c $< -o $@

# Debug compilation rule
$(OBJDIR)/%-debug.o: $(SRCDIR)/%.cpp
	@echo "Compiling $< (debug)..."
	$(CXX) $(DEBUGFLAGS) $(INCLUDES) -MMD -MP -MF $(DEPDIR)/$*-debug.d -c $< -o $@

# ============================================================================
# Utility Targets
# ============================================================================

# Clean all generated files
.PHONY: clean
clean:
	@echo "Cleaning generated files..."
	rm -rf $(OBJDIR) $(BINDIR) $(DEPDIR)
	rm -f *.dat *.log core.*
	@echo "Clean complete."

# Clean only object files (keep binaries)
.PHONY: clean-obj
clean-obj:
	@echo "Cleaning object files..."
	rm -rf $(OBJDIR) $(DEPDIR)

# Show compiler and library information
.PHONY: info
info:
	@echo "Compiler Configuration:"
	@echo "  CXX:      $(CXX)"
	@echo "  CXXFLAGS: $(CXXFLAGS)"
	@echo "  INCLUDES: $(INCLUDES)"
	@echo "  LIBS:     $(LIBS)"
	@echo ""
	@echo "Source Files:"
	@echo "  Shared:   $(SHARED_SOURCES)"
	@echo "  Abundant: $(ABUND_SOURCES)"
	@echo "  Covering: $(COV_SOURCES)"
	@echo ""
	@echo "Targets:"
	@echo "  $(BINDIR)/abundens"
	@echo "  $(BINDIR)/covdens"

# Check dependencies
.PHONY: deps
deps:
	@echo "Checking dependencies..."
	@which $(CXX) > /dev/null || (echo "ERROR: C++ compiler $(CXX) not found" && exit 1)
	@echo "int main(){}" | $(CXX) -x c++ - -lgmp -o /tmp/test_gmp 2>/dev/null || (echo "ERROR: GMP library not found" && exit 1)
	@echo "int main(){}" | $(CXX) -x c++ - -lmpfr -o /tmp/test_mpfr 2>/dev/null || (echo "ERROR: MPFR library not found" && exit 1)
	@echo "int main(){}" | $(CXX) -x c++ - -fopenmp -o /tmp/test_openmp 2>/dev/null || (echo "ERROR: OpenMP not available" && exit 1)
	@rm -f /tmp/test_gmp /tmp/test_mpfr /tmp/test_openmp
	@echo "All dependencies satisfied."

# Show help
.PHONY: help
help:
	@echo "Number Theory Density Calculator Build System"
	@echo "============================================="
	@echo ""
	@echo "Targets:"
	@echo "  all              Build both programs (default)"
	@echo "  abundant         Build abundant number calculator only"
	@echo "  covering         Build covering number calculator only"
	@echo "  debug            Build debug versions of both programs"
	@echo "  debug-abundant   Build debug abundant calculator"
	@echo "  debug-covering   Build debug covering calculator"
	@echo ""
	@echo "Utilities:"
	@echo "  clean            Remove all generated files"
	@echo "  clean-obj        Remove object files only"
	@echo "  deps             Check build dependencies"
	@echo "  info             Show build configuration"
	@echo "  help             Show this help message"
	@echo ""
	@echo "Requirements:"
	@echo "  - C++17 compatible compiler (g++ 7+ recommended)"
	@echo "  - GMP library (GNU Multiple Precision Arithmetic)"
	@echo "  - MPFR library (Multiple Precision Floating-Point)"
	@echo "  - OpenMP support for parallel computation"
	@echo ""
	@echo "Usage Examples:"
	@echo "  make                    # Build both programs"
	@echo "  make abundant           # Build only abundant calculator"
	@echo "  make debug-abundant     # Build debug version"
	@echo "  make clean && make all  # Clean rebuild"

# ============================================================================
# Installation Target (Optional)
# ============================================================================

PREFIX ?= /usr/local
INSTALL_DIR = $(PREFIX)/bin

.PHONY: install
install: all
	@echo "Installing to $(INSTALL_DIR)..."
	install -d $(INSTALL_DIR)
	install -m 755 $(BINDIR)/abundens $(INSTALL_DIR)/
	install -m 755 $(BINDIR)/covdens $(INSTALL_DIR)/
	@echo "Installation complete."

.PHONY: uninstall
uninstall:
	@echo "Uninstalling from $(INSTALL_DIR)..."
	rm -f $(INSTALL_DIR)/abundens $(INSTALL_DIR)/covdens
	@echo "Uninstall complete."

# ============================================================================
# Testing Targets
# ============================================================================

.PHONY: test test-abundant test-covering
test: test-abundant test-covering

test-abundant: $(BINDIR)/abundens
	@echo "Running basic test of abundant number calculator..."
	@echo "Testing with small parameters (zpow=12, numPrimes=100, default precision)..."
	$(BINDIR)/abundens 12 100
	@echo "Abundant number test completed."

test-covering: $(BINDIR)/covdens
	@echo "Running basic test of covering number calculator..."
	@echo "Testing with small parameters (zpow=12, numPrimes=100, default precision)..."
	$(BINDIR)/covdens 12 100
	@echo "Covering number test completed."

# Quick syntax check without full compilation
.PHONY: syntax-check
syntax-check:
	@echo "Checking syntax of all source files..."
	@for src in $(ALL_SOURCES); do \
		echo "  Checking $$src..."; \
		$(CXX) $(CXXFLAGS) $(INCLUDES) -fsyntax-only $$src || exit 1; \
	done
	@echo "Syntax check passed."

# ============================================================================
# Performance and Profiling Targets
# ============================================================================

# Build with profiling enabled
.PHONY: profile
profile: CXXFLAGS += -pg -fno-omit-frame-pointer
profile: $(BINDIR)/abundens-profile $(BINDIR)/covdens-profile

$(BINDIR)/abundens-profile: $(SHARED_OBJECTS) $(ABUND_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

$(BINDIR)/covdens-profile: $(SHARED_OBJECTS) $(COV_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Build optimized version with additional flags
.PHONY: optimized
optimized: CXXFLAGS += -Ofast -funroll-loops -ffast-math
optimized: clean all

# ============================================================================
# Documentation Target
# ============================================================================

.PHONY: doc
doc:
	@echo "Generating documentation..."
	@if command -v doxygen > /dev/null; then \
		doxygen Doxyfile 2>/dev/null || echo "Doxyfile not found, skipping documentation"; \
	else \
		echo "Doxygen not installed, skipping documentation generation"; \
	fi

# ============================================================================
# Include Dependencies
# ============================================================================

# Include dependency files if they exist
-include $(DEPS)

# ============================================================================
# Special Targets
# ============================================================================

# Prevent make from deleting intermediate files
.PRECIOUS: $(OBJDIR)/%.o $(OBJDIR)/%-debug.o

# Declare phony targets
.PHONY: all abundant covering debug debug-abundant debug-covering
.PHONY: clean clean-obj info deps help install uninstall
.PHONY: test test-abundant test-covering syntax-check
.PHONY: profile optimized doc
