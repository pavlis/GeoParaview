/* This object is created from a sta xml file obtained form iris web service.
   */

class DMCWebSta
{
    public:
        DMCWebSta(string fname);
        ~DMCWebSta();
        /* IRIS web service groups stations by net code.
           This method returns the complete list of net codes. */
        list<string>  nets(); 
        /* This returns a vector of station locations for a given net
           code.   */
        vector<SeismicStationLocation> stations(string netcode);
    private:
        mxml_node_t *tree;
};

DMCWebSta::DMCWebSta(string fname)
{
    const string base_error("DMCWebSta constructor:  ");
    FILE *fp;
    fp=fopne(fname.c_str(),"r");
    if(fp==NULL) throw SeisppError(base_error
            + "fopen failure for xml file "+fname);
    tree - mxmlLoadFile(NULL,fp,MXML_NO_CALLBACK);
    flose(fp);
}
DMCWebSta::~DMCWebSta()
{
    mxmlDelete(tree);
}
lixt<string> DMCWebSta::nets()
{
    // This is rough from README file for mxml.   Not sure if the net 
    // is actually the top of the hierarchy
    list<string> result;
    mxml_node_t *node;
    for(node=mxmlFindElement(node,tree,"Network",NULL,NULL,MXML_DESCEND);
        node != NULL; 
        node=mxmlFindELement(node,tree,"Network",NULL,NULL,MXML_DESCEND))
    {
        int whitespacevalue;
        const char *netname = mxmlGetText(node,&whitespacevalue);
        //DEBUG
        cout << "nets method found network="<<node<<endl;
        result.push_back(string(node));
    }
    return(result);
}

