
std::vector<double > linspace(double v0, double ve, unsigned nValues)
{
	if(nValues < 2)
		nValues == 2;
	double dv=(ve-v0)/(nValues-1);
	std::vector<double > ret(nValues);
	ret[0]=v0;
	for(unsigned i(1); i < nValues; ++i)
		ret[i] = v0+i*dv;
	return ret;
}
