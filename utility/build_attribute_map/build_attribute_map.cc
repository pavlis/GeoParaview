#include <iostream>
#include <fstream>
#include <stdio.h>
#include <map>
#include <string>
using namespace std;
string schema_file(string schema_name)
{
	char *ant;
	string antenv("ANTELOPE");
	ant = getenv(antenv.c_str());
	if(ant==NULL)
	{
		cerr << "ANTELOPE environment variable must be defined"<<endl;
		exit(-1);
	}
	return(string(ant)+string("/data/schemas/")+schema_name);
}
void usage()
{
	cerr << "Usage:  build_attribute_map schema_name" <<endl;
	exit(-1);
}	
int main(int argc, char **argv)
{
	string shema_path;
	if(argc!=2)
		usage();
	string schema_name(argv[1]);
	string schema_full_path=schema_file(schema_name);
	ifstream din(schema_full_path.c_str(), ios::in);

	// Skip to first line containing the word "Relation".
	// Assume the file is split so Attribute and Relations 
	// sections are distinct.  Assume Attribute is first with
	// Relation second.  With this assumption skip to the first
	// Attribute line and then work toward the first Relation
	// line.  Parse each Attribute to build amap of type
	// fields.
	char line[1024];
	char *reltoken="Relation";
	char *attkey="Attribute";
	do {
		din.getline(line,1024);
		if(din.eof())
		{
			cerr << "EOF encountered before any Attribute defined for schema "
				<< schema_name << endl
				<< "Probably not a Datascope schema file"<<endl;
			exit(-1);
		}
	} while (strstr(line,attkey)==NULL);

	char creln[40];
	char ctest[20];
	string attribute_name;
	map<string,string> dtypes;
	do {
		// Just look for a fixed set of keywords and build the map
		if(strstr(line,attkey)!=NULL)
		{
			sscanf(line,"%s%s",ctest,creln);
			attribute_name=string(creln);
		}
		else if(strstr(line,"String")!=NULL)
			dtypes[attribute_name]="string";
		else if(strstr(line,"Real")!=NULL)
			dtypes[attribute_name]="real";
		else if(strstr(line,"Integer")!=NULL)
			dtypes[attribute_name]="int";
		else if(strstr(line,"Time")!=NULL)
			dtypes[attribute_name]="real";
		din.getline(line,1024);

	} while (!din.eof());
	// rewind and skip to the first Relation line
	//din.seekg(0);
	// could not make rewind work, close and open
	din.close();
	din.open(schema_full_path.c_str(), ios::in);
	do {
		din.getline(line,1024);
	}
	while (strstr(line,reltoken)==NULL);

// Assume all the Relation definitions are together.  This will fail
// if they are mixed with Attribute definitions
	cout << schema_name << " Tbl&{"<<endl;
	do {
		sscanf(line,"%s%s",ctest,creln);
		if(strcmp(ctest,reltoken))
		{
			cerr << "Fatal parsing error"<<endl
			   << "Expected this line ->"<<line<<endl
			   << "To begin with the word->Relation<-"<<endl;
			exit(-1);
		}
		string table(creln);
		do {
			din.getline(line,1024);
		} while ((strstr(line,"Fields")==NULL) && !din.eof());
		if(din.eof()) 
		{
			cerr<<"Warning premature exit searching for Fields keyword"
				<< "Last relation definitions may have been dropped"
				<< endl;
			exit(0);
		}
		// line should now contain the Fields ( ) component
		
		string sline(line);
		const string white(" \t\n)"); // Need ) for cases like lddate)
		int current, end_current;
		int end_line;
		end_line = sline.find(")");
		if(end_line<0)
		{
		// Have to hunt for an end in this condition as the ) can
		// be on a following line
		    string temp;
		    do
		    {
			din.getline(line,1024);
			if(din.eof())
			{
			   cerr << "Parse error:  Hit EOF search for ) to close Fields line"
				<< endl;
			  exit(-1);
			}
			temp = string(line);
			sline += temp;
			end_line = sline.find(")");
		    } while(end_line<0);
		}
		current = sline.find("(");
		if(current<0)
		{
			cerr << "Parse error:  cannot find ( in Fields line"
				<< endl;
			exit(-1);
		}
		++current;
		do {
			current=sline.find_first_not_of(white,current);
			if(current<0) break; // needed to break on )
			end_current=sline.find_first_of(white,current);
			if(end_current==end_line)
				attribute_name.assign(sline,current,
					end_current-current-1);
			else
				attribute_name.assign(sline,current,
					end_current-current);
			if(attribute_name == ")") 
				break;
			else
			{
				cout << table <<"."<<attribute_name<< " "
				<< attribute_name << " "
				<< table << " "
				<< dtypes[attribute_name] << endl;
			}
			current=end_current;
		}
		while(end_current<end_line);
		// Hunt for the next Relation line
		do {
			din.getline(line,1024);
			if(din.eof()) 
			{
				// this is the only exit -- ugly
				cout << "}" << endl;
				exit(0);
			}
		} while (strstr(line,reltoken)==NULL );

	} while (1);
}
