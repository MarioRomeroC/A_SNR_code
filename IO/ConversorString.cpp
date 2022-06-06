/**
    Some functions involving strings. Added as 'mandatory' utility and to avoid rewritting again and again.

    This is an adapted version of part the code used for analysis in Romero & Ascasibar 2018.

    Mario Romero Spring 2015
**/
#ifndef CONVERSOR_STRING_CPP
#define CONVERSOR_STRING_CPP

#include <sstream> //Only used here

string tostring(int a){
	//int->string conversion.
	stringstream conv;
	conv<<a;
	return conv.str();
}
string tostring(my_float a){
	//float->string conversion.
	stringstream conv;
	conv<<a;
	return conv.str();
}
my_float todouble(string A){
	//string->float conversion
	my_float n;
	stringstream conv;
	conv<<A;
	conv>>n;
	return n;
}
my_float todouble(char *A[]){
    //string->float conversion
	my_float n;
	stringstream conv;
	conv<<A;
	conv>>n;
	return n;
}
bool identify(string word,string phrase){
    //Search a word in a string
	if(phrase.find(word) != string::npos){ //It does.
		return true;}
	return false;
}

string extract(string phrase, string start, string end){
    //Extract a string between start (WITHOUT COUNTING IT) and end.
	size_t startpos = phrase.find(start) + start.length(); //Issue here is that position is given alonside the string 'start', and so you have to include the addition
	size_t endpos = phrase.find(end);
	size_t length = (endpos-startpos);

	return phrase.substr(startpos,length);

}

#endif
