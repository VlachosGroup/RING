#ifndef CLONABLE_H
#define CLONABLE_H

//#include <iostream>
//#include <fstream>
//#include <string>
//#include <cstring>
//#include <sstream>
//#include <vector>
//#include <list>
//#include <set>
//#include <map>
//#include <queue>
//#include <deque>
//#include <algorithm>
//using namespace std;

class clonable {

	public:
		virtual ~clonable() {}
		virtual clonable* clone() = 0;
};

#endif
