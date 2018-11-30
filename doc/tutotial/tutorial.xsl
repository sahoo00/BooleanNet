<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  version="1.0">
  <xsl:import href="/usr/share/docbook-xsl/htmlhelp/htmlhelp.xsl"/>
  <xsl:param name="generate.legalnotice.link" select="1"/>
  <xsl:param name="suppress.navigation" select="0"/>
  <xsl:param name="admon.graphics" select="1"/>
  <xsl:param name="admon.graphics.path">images/</xsl:param>
  <xsl:param name="htmlhelp.chm" select="'tutorial.chm'"/>
  <xsl:param name="htmlhelp.hhc.binary" select="0"/>
  <xsl:param name="htmlhelp.hhc.folders.instead.books" select="0"/>
  <xsl:param name="toc.section.depth" select="4"/>
  <!--
  <xsl:param name="html.stylesheet" select="'tutorial.css'"/>
  <xsl:template name="user.header.navigation">
    <hr></hr>
    <p>StepMiner Documentation</p>
    <hr></hr>
  </xsl:template>

  <xsl:template name="user.footer.navigation">
    <hr></hr>
    <p>StepMiner: Version 1.0 Alpha</p>
    <hr></hr>
  </xsl:template>
  -->
</xsl:stylesheet>
