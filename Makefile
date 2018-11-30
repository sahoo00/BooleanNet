release: clean release1 release2 macosx archive
release1:
	mkdir -p StepMiner-1.0
	cp dist/lib/tools.jar StepMiner-1.0/stepminer-1.0.jar
	cp demoexamples/yeast.pcl StepMiner-1.0
	cp demoexamples/serum.pcl StepMiner-1.0
	cp demoexamples/hongjuan-10uM.pcl StepMiner-1.0
	cp gene_* StepMiner-1.0
	cp -r images StepMiner-1.0
	cp -r doc StepMiner-1.0
	echo java -Xms64m -Xmx512m -jar stepminer-1.0.jar --gui > StepMiner-1.0/stepminer.bat
	echo java -Xms64m -Xmx512m -jar stepminer-1.0.jar --gui > StepMiner-1.0/stepminer
	echo "#!/bin/sh" > StepMiner-1.0/stepminer.sh
	echo java -Xms64m -Xmx512m -jar stepminer-1.0.jar --gui >> StepMiner-1.0/stepminer.sh
release2:
	mkdir -p StepMiner-1.0-java-1.4
	cp build1.4/stepminer-1.0.jar StepMiner-1.0-java-1.4/stepminer-1.0.jar
	cp demoexamples/yeast.pcl StepMiner-1.0-java-1.4
	cp demoexamples/serum.pcl StepMiner-1.0-java-1.4
	cp demoexamples/hongjuan-10uM.pcl StepMiner-1.0-java-1.4
	cp gene_* StepMiner-1.0-java-1.4
	cp -r images StepMiner-1.0-java-1.4
	cp -r doc StepMiner-1.0-java-1.4
	echo java -Xms64m -Xmx512m -jar stepminer-1.0.jar --gui > StepMiner-1.0-java-1.4/stepminer.bat
	echo java -Xms64m -Xmx512m -jar stepminer-1.0.jar --gui > StepMiner-1.0-java-1.4/stepminer
macosx:
	cp -r StepMiner-1.0-java-1.4 StepMiner-1.0-osx
	cp -r dist/StepMiner.app StepMiner-1.0-osx
	cp build1.4/stepminer-1.0.jar StepMiner-1.0-osx/StepMiner.app/Contents/Resources/Java/tools.jar
archive:
	tar cvzf StepMiner-1.0.tar.gz StepMiner-1.0
	tar cvzf StepMiner-1.0-java-1.4.tar.gz StepMiner-1.0-java-1.4
	zip -r StepMiner-1.0.zip StepMiner-1.0
	zip -r StepMiner-1.0-java-1.4.zip StepMiner-1.0-java-1.4
	zip -r StepMiner-1.0-osx.zip StepMiner-1.0-osx
zip:
	zip -r StepMiner.zip doc images build.xml LICENSE INSTALL README \
	dist tools stepminer stepminer.bat \
	manifest setupmac demoexamples gene_association.goa_human \
	gene_association.sgd gene_ontology.obo
ozip:
	zip -r StepMiner.zip doc images build.xml LICENSE INSTALL README \
	dist tools stepminer stepminer.bat Makefile \
	manifest setupmac demoexamples gene_association.goa_human \
	gene_association.sgd gene_ontology.obo build1.4 lib
clean:
	rm -fr StepMiner-1.0
	rm -fr StepMiner-1.0-java-1.4
	rm -fr StepMiner-1.0.zip
	rm -fr StepMiner-1.0.tar.gz
	rm -fr StepMiner-1.0-java-1.4.zip
	rm -fr StepMiner-1.0-java-1.4.tar.gz
	rm -fr StepMiner-1.0-osx
	rm -fr StepMiner-1.0-osx.zip
	rm -fr StepMiner.zip
