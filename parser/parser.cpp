#include "parser.h"

static	void	add_x_id(vector<vector<VALUE_TYPE>> &slae)
{
	vector<VALUE_TYPE> x;

	for (int i = 1; i <= slae.size(); i++)
		x.push_back(i);

	// id = 0 relates to b column
	x.push_back(0);

	slae.push_back(x);
}

static void 	remove_outer_spaces(string *line)
{
	int pos1 = line->find_first_not_of(' ');
	line->erase(0, pos1);

	int pos2 = line->find_last_not_of(' ');
	line->erase(pos2 + 1, line->length());
}

vector<vector<VALUE_TYPE>> 	parsing_file()
{
	vector<vector<VALUE_TYPE>> slae;

	ifstream	slae_file(FILE_NAME);
	if (!slae_file.is_open())
	{
		cout << endl << RED << "[WARN] File " << FILE_NAME << " not open!" << RESET << endl;
		exit(0);
	}

	string	line;
	while (getline(slae_file, line))
	{
		vector<VALUE_TYPE> tmp;

		if (line.find_first_of("0123456789") == string::npos || line.find_first_of("*") != string::npos)
			continue;

		remove_outer_spaces(&line);

		int pos;
		while (line.find_first_of("0123456789") != string::npos)
		{
			tmp.push_back(typeid(VALUE_TYPE) == typeid(double) ? stod(line) : stof(line));
			if ((pos = line.find(' ')) == string::npos)
				break;
			line.erase(0, pos + 1);
			line.erase(0, line.find_first_not_of(' '));
		}

		slae.push_back(tmp);
	}

	/* удаляем в итерационных методах хранение номера столбцов */
//	add_x_id(slae);

	slae_file.close();

	return slae;
}

