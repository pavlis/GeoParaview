#include <string>
#include <ostream>
using namespace std;
int main(int argc, char **argv)
{
	// tested string stmp="test1 test2	test3	test4";
	// tested string stmp="  test1 test2	test3	test4";
	string stmp="test1 test2	test3	test4  ";
	// tested string stmp="  test1 test2	test3	test4  ";
	int current, end_current;
	string db_attribute_name,db_table_name,internal_name,s4;
	const string white(" \t\n");

        if((current=stmp.find_first_not_of(white,0)) != 0)
        {
                stmp.erase(0,current);
                current = 0;
        }
        end_current=stmp.find_first_of(white,current);
        db_attribute_name.assign(stmp,current,end_current-current);
cout << "token 1 ="<<db_attribute_name<<endl;
        current = stmp.find_first_not_of(white,end_current);
        end_current=stmp.find_first_of(white,current);
        db_table_name.assign(stmp,current,end_current-current);
cout << "token 2 ="<<db_table_name<<endl;
        current = stmp.find_first_not_of(white,end_current);
        end_current=stmp.find_first_of(white,current);
        internal_name.assign(stmp,current,end_current-current);
cout << "token 3 ="<<internal_name<<endl;
        current = stmp.find_first_not_of(white,end_current);
        end_current=stmp.find_first_of(white,current);
	if(end_current<0)end_current=stmp.length();
        s4.assign(stmp,current,end_current-current);
cout << "token 4 ="<<s4<<endl;
}

	

	
