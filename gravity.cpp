#include <cstdio>
#include <cmath>

// Units are in kilograms and meters. ---------------------------------------

// Vector functions ---------------------------------------------------------

void vec3Set(double x, double y, double z, double r[3])
{
    r[0] = x;
    r[1] = y;
    r[2] = z;
}

void vec3Divide(double v0[3], double d, double r[3])
{
    for(int m = 0; m < 3; m++) {
        r[m] = v0[m] / d;
    }
}

void vec3Multiply(double v0[3], double d, double r[3])
{
    for(int m = 0; m < 3; m++) {
        r[m] = v0[m] * d;
    }
}

void vec3Add(double v0[3], double v1[3], double r[3])
{
    for(int m = 0; m < 3; m++) {
        r[m] = v0[m] + v1[m];
    }
}

void vec3Subtract(double v0[3], double v1[3], double r[3])
{
    for(int m = 0; m < 3; m++) {
        r[m] = v0[m] - v1[m];
    }
}

double vec3Length(double v[3])
{
    return sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}


// Physics functions --------------------------------------------------------

const double gravitationalConstant = 6.673 * 1e-11; // 10-11 N*m2/kg2;

double gravity(double mass, double distance)
{
    return gravitationalConstant * mass / (distance * distance);
}


// Simulation parameters and functions --------------------------------------

const int res = 1000;

// Mass of moon rock in kg/(m^3)
const double moonMass = 3.346 * 100 * 100 * 100 / 1000; // 3.346 g/cm3
const double moonRadius = 1737.1 * 1000; // 1737.1 km

// Convert unit cube coordinate to real-world meters
// In this case, a little larger than the diameter of the moon
double unitConvert(double x)
{
    // return x * moonRadius * 1.1;
    return x * 2 * moonRadius * 1.1;
}

// Mass in volume
// In this case, moon mass in kg/(m^3)
double massInVolume(double w, double h, double d)
{
    return moonMass * w * h * d;
}

// Moon
bool sampleWithinShape(double sample[3])
{
    double distanceToCenter = sqrtf(sample[0] * sample[0] + sample[1] * sample[1] + sample[2] * sample[2]);

    return distanceToCenter < moonRadius;
}


// Simulation loop ----------------------------------------------------------

int main()
{
    double d = unitConvert(2.0 / res);
    double mass = massInVolume(d, d, d);

    double myPosition[3];
    // Standing on the surface
    vec3Set(0.0, moonRadius, 0.0, myPosition);
    double myMass = 106.594;

    double finalVector[3];
    vec3Set(0.0, 0.0, 0.0, finalVector);

    for(int k = 0; k < res; k++) {

        double accumulatedSlice[3];
        vec3Set(0.0, 0.0, 0.0, accumulatedSlice);

        for(int j = 0; j < res; j++) {

            double accumulatedRow[3];
            vec3Set(0.0, 0.0, 0.0, accumulatedRow);

            for(int i = 0; i < res; i++) {

                double sample[3];
                sample[0] = unitConvert(2.0 / res * (i - res/2.0 + .5));
                sample[1] = unitConvert(2.0 / res * (j - res/2.0 + .5));
                sample[2] = unitConvert(2.0 / res * (k - res/2.0 + .5));

                if(sampleWithinShape(sample)) {

                    double difference[3];
                    vec3Subtract(myPosition, sample, difference);

                    double distance = vec3Length(difference);

                    double direction[3];
                    vec3Divide(difference, distance, direction);

                    double g = gravity(mass, distance);
                    double acceleration[3];
                    vec3Multiply(direction, g, acceleration);

                    vec3Add(accumulatedRow, acceleration, accumulatedRow);
                }
            }

            vec3Add(accumulatedSlice, accumulatedRow, accumulatedSlice);
        }

        vec3Add(finalVector, accumulatedSlice, finalVector);
    }

    printf("acceleration is %g m*s^-2\n", vec3Length(finalVector));
}
