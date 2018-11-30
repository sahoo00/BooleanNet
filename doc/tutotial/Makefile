all: tutorial.xml
	xsltproc --nonet tutorial.xsl tutorial.xml
	#xsltproc --nonet /usr/share/docbook-xsl/htmlhelp/htmlhelp.xsl tutorial.xml
chm: htmlhelp.hhp
	"C:/Program Files/HTML Help Workshop/hhc.exe" htmlhelp.hhp 
