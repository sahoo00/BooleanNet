/*

Copyright (c) 2006, the Board of Trustees of Leland
Stanford Junior University.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the
distribution.

    * Neither the name of Stanford University nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/*
 Author: Debashis Sahoo <sahoo@stanford.edu>
 Author: David Dill <dill@cs.stanford.edu>
 */

package tools.parsers.Expression;

import java_cup.runtime.*;
      
%%
   
/* Generated class will be Lexer */
%class Lexer

/* Keep line and column info */
%line
%column

%cup
   
/*
  Declarations
   
  Code between %{ and %}, both of which must be at the beginning of a
  line, will be copied letter to letter into the lexer class source.
  Here you declare member variables and functions that are used inside
  scanner actions.  
*/

%{   
    StringBuffer string = new StringBuffer();

    /* To create a new java_cup.runtime.Symbol with information about
       the current token, the token will have no value in this
       case. */
    private Symbol symbol(int type) {
        return new Symbol(type, yyline, yycolumn);
    }
    
    /* Also creates a new java_cup.runtime.Symbol with information
       about the current token, but this object has a value. */
    private Symbol symbol(int type, Object value) {
        return new Symbol(type, yyline, yycolumn, value);
    }
%}

LineTerminator = \r|\n|\r\n
   
WhiteSpace     = {LineTerminator} | [ \t\f]

/* This is a gene, probeset ID, virtual gene name, function name, etc. */
ID = ([A-Za-z0-9_]+)

/* Decimal integer or floating point.  Can be 012 or 1. or 1.1 or .03 but not "." */

Number = (([0-9]+("."[0-9]*)?)|([0-9]*"."[0-9]+))

%state STR

%%
   
<YYINITIAL> {
   
    ";"                { return symbol(sym.SEMI); }
    ","                { return symbol(sym.COMMA); }
    "=="                { return symbol(sym.EQ); }   
    "="                { return symbol(sym.ASSN); }   
    "<"                { return symbol(sym.LT); }   
    ">"                { return symbol(sym.GT); }   
    "<="                { return symbol(sym.LEQ); }   
    ">="                { return symbol(sym.GEQ); }   
    "+"                { return symbol(sym.PLUS); }
    "-"                { return symbol(sym.MINUS); }
    "*"                { return symbol(sym.TIMES); }
    "/"                { return symbol(sym.DIVIDE); }
    "&"                { return symbol(sym.AND); }
    "|"                { return symbol(sym.OR); }
    "->"                { return symbol(sym.IMPLIES); }
    "!"                { return symbol(sym.NOT); }
    "("                { return symbol(sym.LPAREN); }
    ")"                { return symbol(sym.RPAREN); }

    /* Builtin keywords.  These are recognized before ID */

    "phenotype" { return symbol(sym.PHENOTYPE); }
    "search" { return symbol(sym.SEARCH); }
    "csearch" { return symbol(sym.CSEARCH); }
    "rsearch" { return symbol(sym.RSEARCH); }
    "threshold" { return symbol(sym.THRESHOLD); }
    "upper" { return symbol(sym.UPPER); }
    "lower" { return symbol(sym.LOWER); }
    "high" { return symbol(sym.HIGH); }
    "low" { return symbol(sym.LOW); }
    "med" { return symbol(sym.MED); }
    "virtual" { return symbol(sym.VIRTUAL); }

    {Number}  { return symbol(sym.NUMBER, new Double(yytext()));}
   
    {ID}      { return symbol(sym.ID, new String(yytext())); }
    [$]{ID}      { return symbol(sym.TMPID, new String(yytext())); }
   
    {WhiteSpace}       { /* just skip what was found, do nothing */ }   
    \"                 { string.setLength(0); yybegin(STR); }
}

<STR> {
  \"                             { yybegin(YYINITIAL); 
                                   return symbol(sym.STRING, 
                                   string.toString()); }
  [^\n\r\"\\]+                   { string.append( yytext() ); }
  \\t                            { string.append('\t'); }
  \\n                            { string.append('\n'); }

  \\r                            { string.append('\r'); }
  \\\"                           { string.append('\"'); }
  \\                             { string.append('\\'); }
}

/* No token was found for the input so through an error.  Print out an
   Illegal character message with the illegal character that was found. */
[^]                    { throw new Error("Illegal character <"+yytext()+">"); }
