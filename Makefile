
JAVAC = javac
JAVAC_FLAGS = -g -classpath .:hb15.zip

PACKAGES = \
	common \
	alignCompress \
	fuzzyLZ

JAVA_DIRS       := $(subst .,/,$(PACKAGES))
JAVA_SRC        := $(foreach dir, $(JAVA_DIRS), $(wildcard $(dir)/*.java))

JAVA_OBJS       := $(JAVA_SRC:.java=.class)



all: $(JAVA_OBJS)



%.class : %.java
	${JAVAC} $(JAVAC_FLAGS) $<

clean:
	find . \( -name '*~' -o -name '*.class' \) -print | xargs rm
