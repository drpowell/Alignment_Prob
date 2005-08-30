
JAVAC = javac
#JAVAC_FLAGS = -g
JAVAC_FLAGS = -O -target 1.2 -source 1.2

PACKAGES = \
	common \
	alignCompress \
	fuzzyLZ

COMMON_DIRS     = common
COMMON_SRC      := $(foreach dir, $(COMMON_DIRS), $(wildcard $(dir)/*.java))
COMMON_OBJS     := $(COMMON_SRC:.java=.class)

FUZZYLZ_SRC  := $(wildcard fuzzyLZ/*.java)
FUZZYLZ_OBJS := $(FUZZYLZ_SRC:.java=.class)

ALIGNCOMPRESS_SRC  := $(wildcard alignCompress/*.java)
ALIGNCOMPRESS_OBJS := $(ALIGNCOMPRESS_SRC:.java=.class)


all: fuzzyLZ.jar alignCompress.jar

fuzzyLZ.jar: $(COMMON_OBJS) $(FUZZYLZ_OBJS)
	jar cmf fuzzyLZ/myManifest fuzzyLZ.jar common/*.class fuzzyLZ/*.class

fuzzyLZ-src.tar.gz: $(COMMON_OBJS) $(FUZZY_OBJS) fuzzyLZ.jar COPYRIGHT README.FuzzyLZ Makefile.fuzzyLZ fuzzyLZ/myManifest smooth.pl fuzzyLZ-params.txt
	./tar.pl fuzzyLZ-src.tar.gz fuzzyLZ fuzzyLZ.jar $(COMMON_SRC) $(FUZZYLZ_SRC) COPYRIGHT README.FuzzyLZ Makefile.fuzzyLZ=Makefile fuzzyLZ/myManifest smooth.pl fuzzyLZ-params.txt

alignCompress.jar: $(COMMON_OBJS) $(ALIGNCOMPRESS_OBJS)
	jar cmf alignCompress/myManifest alignCompress.jar common/*.class alignCompress/*.class

alignCompress-src.tar.gz: alignCompress.jar $(COMMON_SRC) $(ALIGNCOMPRESS_SRC) COPYRIGHT README.alignCompress
	./tar.pl alignCompress-src.tar.gz alignCompress alignCompress.jar $(COMMON_SRC) $(ALIGNCOMPRESS_SRC) COPYRIGHT README.alignCompress

%.class : %.java
	${JAVAC} $(JAVAC_FLAGS) $<

clean:
	find . \( -name '*~' -o -name '*.class' \) -print | xargs rm
	rm -f fuzzyLZ.jar alignCompress.jar
