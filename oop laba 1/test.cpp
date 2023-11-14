#include "catch.hpp"
#include "function.h"

using namespace std;
TEST_CASE("standard distribution")
{
	HuberDistribution* HD = new HuberDistribution();
	HD->v = 0.5;
	HD->K = K(HD->v);
	HD->scale = 1.;
	HD->shift = 0.;
	REQUIRE(Dksi_huber(HD) == Approx(8.08).epsilon(0.01));
	REQUIRE(kurtosis_huber(HD) == Approx(2.94).epsilon(0.01));
	REQUIRE(P(HD) == Approx(0.214).epsilon(0.01));
	REQUIRE(Huber(0., HD) == Approx(0.223).epsilon(0.01));
	HD->v = 0.75;
	HD->K = K(HD->v);
	REQUIRE(Dksi_huber(HD) == Approx(3.71).epsilon(0.01));
	REQUIRE(kurtosis_huber(HD) == Approx(2.75).epsilon(0.01));
	REQUIRE(P(HD) == Approx(0.405).epsilon(0.01));
	REQUIRE(Huber(0., HD) == Approx(0.296).epsilon(0.01));
	HD->v = 1;
	HD->K = K(HD->v);
	REQUIRE(Dksi_huber(HD) == Approx(2.24).epsilon(0.01));
	REQUIRE(kurtosis_huber(HD) == Approx(2.37).epsilon(0.01));
	REQUIRE(P(HD) == Approx(0.585).epsilon(0.01));
	REQUIRE(Huber(0., HD) == Approx(0.342).epsilon(0.01));
	HD->v = 1.5;
	HD->K = K(HD->v);
	REQUIRE(Dksi_huber(HD) == Approx(1.31).epsilon(0.01));
	REQUIRE(kurtosis_huber(HD) == Approx(1.30).epsilon(0.01));
	REQUIRE(P(HD) == Approx(0.834).epsilon(0.01));
	REQUIRE(Huber(0., HD) == Approx(0.384).epsilon(0.01));
	HD->v = 2;
	HD->K = K(HD->v);
	REQUIRE(Dksi_huber(HD) == Approx(1.08).epsilon(0.01));
	REQUIRE(kurtosis_huber(HD) == Approx(0.51).epsilon(0.01));
	REQUIRE(P(HD) == Approx(0.946).epsilon(0.01));
	REQUIRE(Huber(0., HD) == Approx(0.396).epsilon(0.01));
	HD->v = 2.5;
	HD->K = K(HD->v);
	REQUIRE(Dksi_huber(HD) == Approx(1.02).epsilon(0.01));
	REQUIRE(kurtosis_huber(HD) == Approx(0.16).epsilon(0.1));
	REQUIRE(P(HD) == Approx(0.986).epsilon(0.01));
	REQUIRE(Huber(0., HD) == Approx(0.398).epsilon(0.01));
	HD->v = 3;
	HD->K = K(HD->v);
	REQUIRE(Dksi_huber(HD) == Approx(1.00).epsilon(0.01));
	REQUIRE(kurtosis_huber(HD) == Approx(0.04).epsilon(0.01));
	REQUIRE(P(HD) == Approx(0.997).epsilon(0.01));
	REQUIRE(Huber(0., HD) == Approx(0.399).epsilon(0.01));
}
