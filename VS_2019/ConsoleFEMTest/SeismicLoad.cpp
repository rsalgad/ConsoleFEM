#include "pch.h"

SeismicLoad::SeismicLoad()
{
}

std::vector<double> SeismicLoad::GetTime()
{
	return _time;
}

void SeismicLoad::SetTime(std::vector<double> time)
{
	_time = time;
}

std::vector<std::vector<double>> SeismicLoad::GetRecords()
{
	return _records;
}

std::vector<char> SeismicLoad::GetDirections()
{
	return _directions;
}

void SeismicLoad::SetRecordX(std::vector<double> record)
{
	_records.push_back(record);
	_directions.push_back('x');
}

void SeismicLoad::SetRecordY(std::vector<double> record)
{
	_records.push_back(record);
	_directions.push_back('y');
}

void SeismicLoad::SetRecordZ(std::vector<double> record)
{
	_records.push_back(record);
	_directions.push_back('z');
}

double SeismicLoad::LoadFromTime(double t, char dir) const{
	int index = GetIndexOfDirection(dir);
	double g = 9.8 * 1000;
	if (index != -1) {
		int nextIndex = 0;
		int prevIndex = 0;
		for (int i = 0; i < _time.size(); i++) {
			if (t <= _time[i]) {
				nextIndex = i;
				if (nextIndex != 0) {
					prevIndex = i - 1;
				}
				break;
			}
			else {
				nextIndex = -1; // time is beyond the timespan of the seismic history
			}
		}

		if (nextIndex != -1) {
			double nextLoad = _records[index][nextIndex] * g;
			double nextTime = _time[nextIndex];
			double prevLoad, prevTime;
			if (nextIndex == 0) {
				prevLoad = 0;
				prevTime = 0;
			}
			else {
				prevLoad = _records[index][prevIndex] * g;
				prevTime = _time[prevIndex];
			}

			return (prevLoad * (nextTime - t) + nextLoad * (t - prevTime)) / (nextTime - prevTime);
		}
		else {
			return 0;
		}
	}
	else {
		return 0;
	}

}

int SeismicLoad::GetIndexOfDirection(char dir) const {
	int index = -1;
	for (int j = 0; j < _directions.size(); j++) {
		if (_directions[j] == dir) {
			index = j;
		}
	}
	return index;
}

std::string SeismicLoad::GetType()
{
	return std::string();
}



Matrix SeismicLoad::GetSeismicLoadVector(const SeismicLoad& sLoad, Matrix& FInc, const double* t) {
	Matrix ans(FInc.GetDimX(), 1);

	double xLoad, yLoad, zLoad;
	if (sLoad.GetIndexOfDirection('x') != -1) {
		xLoad = sLoad.LoadFromTime(*t, 'x');
	}
	if (sLoad.GetIndexOfDirection('y') != -1) {
		yLoad = sLoad.LoadFromTime(*t, 'y');
	}
	if (sLoad.GetIndexOfDirection('z') != -1) {
		zLoad = sLoad.LoadFromTime(*t, 'z');
	}

	for (int i = 0; i < FInc.GetDimX(); i++) {
		if (FInc.GetMatrixDouble()[i][0] == 1) {
			ans.GetMatrixDouble()[i][0] = xLoad;
		}
		else if (FInc.GetMatrixDouble()[i][0] == 2)
		{
			ans.GetMatrixDouble()[i][0] = yLoad;
		}
		else if (FInc.GetMatrixDouble()[i][0] == 3) {
			ans.GetMatrixDouble()[i][0] = zLoad;
		}
	}
	return ans;
}

SeismicLoad::~SeismicLoad()
{
}
