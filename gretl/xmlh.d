/*
   Minimal D declarations for the parts of libxml2 that the ported gretl
   code touches directly (gretl_xml.h includes the libxml2 headers). libxml2
   is an external library, so only what is needed is declared, as extern (C)
   types; add more here as ported files need it.

   xmlNode is declared only as far as the "next" member (the leading part of
   the real struct, which is how libxml2's own headers lay it out); code
   must only ever handle xmlNode pointers obtained from libxml2, never
   allocate or copy an xmlNode.
*/

module gretl.xmlh;


alias xmlChar = ubyte;
alias xmlElementType = int;

struct xmlDoc;
alias xmlDocPtr = xmlDoc *;

struct xmlNode {
    void * _private;
    xmlElementType type;
    const(xmlChar) * name;
    xmlNode * children;     /* xmlChildrenNode */
    xmlNode * last;
    xmlNode * parent;
    xmlNode * next;
    /* ... remaining members not declared ... */
}
alias xmlNodePtr = xmlNode *;

/* xmlError, as far as the "line" member */
struct xmlError {
    int domain;
    int code;
    char * message;
    int level;              /* xmlErrorLevel */
    char * file;
    int line;
    /* ... remaining members not declared ... */
}

/* xmlParserOption values used by gretl */
enum XML_PARSE_NOBLANKS = 1 << 8;
enum XML_PARSE_HUGE = 1 << 19;

/* DTD validation (used by gretl_func for checking function packages) */
struct xmlDtd;
alias xmlDtdPtr = xmlDtd *;
alias xmlValidityErrorFunc = extern (C) void function(void * ctx, const(char) * msg, ...);
alias xmlValidityWarningFunc = extern (C) void function(void * ctx, const(char) * msg, ...);

struct xmlValidCtxt {
    void * userData;                   /* user specific data block */
    xmlValidityErrorFunc error;        /* the callback in case of errors */
    xmlValidityWarningFunc warning;    /* the callback in case of warning */
    /* ... remaining members not declared ... */
}
alias xmlValidCtxtPtr = xmlValidCtxt *;

/* libxml2 string helpers and the allocator hook, used with xmlNode data */
extern (C) {
    int xmlStrcmp(const(xmlChar) * str1, const(xmlChar) * str2);
    xmlChar * xmlGetProp(const(xmlNode) * node, const(xmlChar) * name);
    xmlChar * xmlNodeGetContent(const(xmlNode) * cur);
    xmlChar * xmlNodeListGetString(xmlDocPtr doc, const(xmlNode) * list, int inLine);
    xmlChar * xmlNodeListGetRawString(const(xmlDoc) * doc, const(xmlNode) * list, int inLine);
    int xmlIsBlankNode(const(xmlNode) * node);
    xmlDocPtr xmlReadFile(const(char) * filename, const(char) * encoding, int options);
    xmlDocPtr xmlReadMemory(const(char) * buffer, int size, const(char) * url,
                            const(char) * encoding, int options);
    xmlNodePtr xmlDocGetRootElement(const(xmlDoc) * doc);
    void xmlFreeDoc(xmlDocPtr cur);
    xmlDtdPtr xmlParseDTD(const(xmlChar) * ExternalID, const(xmlChar) * SystemID);
    void xmlFreeDtd(xmlDtdPtr cur);
    xmlValidCtxtPtr xmlNewValidCtxt();
    void xmlFreeValidCtxt(xmlValidCtxtPtr cur);
    int xmlValidateDtd(xmlValidCtxtPtr ctxt, xmlDocPtr doc, xmlDtdPtr dtd);
    const(xmlError) * xmlGetLastError();
    void xmlInitParser();
    void xmlCleanupParser();
    version (Windows) {
        /* xmlFree is a data export of libxml2, which a DLL import library does not resolve; its default value is xmlMemFree */
        void xmlMemFree(void *);
    } else {
        extern __gshared void function(void *) xmlFree;
    }
}

version (Windows) alias xmlFree = xmlMemFree;
