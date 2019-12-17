#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>

#include "vec3.h"
#include "time.h"

#define G 4*M_PI*M_PI
#define c  63239.7263

using namespace std;


class CelestialBody
{
public:
    vec3 position;
    vec3 velocity;
    vec3 force;
    double mass;

    CelestialBody(vec3 pos, vec3 vel, double mass_) {
    position = pos;
    velocity = vel;
    mass = mass_;
    }
    CelestialBody(double x, double y, double z, double vx, double vy, double vz, double mass_) {
    position = vec3(x,y,z);
    velocity = vec3(vx,vy,vz);
    mass = mass_;
    }
    void resetForce() {
    force.zeros();
    }
};



class SolarSystem
{
public:
    SolarSystem();
    CelestialBody &createCelestialBody(vec3 position, vec3 velocity, double mass);
    void calculateEnergyMomentum();
    void calculateGravitationalForce();
    void calculateGravitationalForce_relativistic();
    int numberOfBodies() const;

    double totalEnergy() const;
    double potentialEnergy() const;
    double kineticEnergy() const;
    void writeToFile(std::string filename);
    void writeToFile_Energy_Momentum(std::string filename, double, double, double, double);
    vec3 angularMomentum() const;
    std::vector<CelestialBody> &bodies();

private:
    std::vector<CelestialBody> m_bodies;
    vec3 m_angularMomentum;
    std::ofstream m_file;
    std::ofstream n_file;
    double m_kineticEnergy;
    double m_potentialEnergy;
    double L_x, L_y, L_z;
    double L_x_Merc, L_y_Merc, L_z_Merc, L_Merc;
};



SolarSystem::SolarSystem() :
    m_kineticEnergy(0),
    m_potentialEnergy(0)
{
}

CelestialBody& SolarSystem::createCelestialBody(vec3 position, vec3 velocity, double mass) {
    m_bodies.push_back( CelestialBody(position, velocity, mass) );
    return m_bodies.back(); // Return reference to the newest added celstial body
}

void SolarSystem::calculateEnergyMomentum()
{
    m_kineticEnergy = 0;
    m_potentialEnergy = 0;
    m_angularMomentum.zeros();
    L_x = 0, L_y = 0, L_z = 0; //components of angular momentum
    

    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body1 = m_bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++) {
            CelestialBody &body2 = m_bodies[j];
            vec3 deltaRVector = body1.position - body2.position;
            double dr = deltaRVector.length();
            // Calculate potential energy here
            m_potentialEnergy += -G*body1.mass*body2.mass/dr;
           
        }
         L_x += body1.mass*(body1.position.y()*body1.velocity.z() - body1.position.z()*body1.velocity.y());
         L_y += body1.mass*(body1.position.z()*body1.velocity.x() - body1.position.x()*body1.velocity.z());
         L_z += body1.mass*(body1.position.x()*body1.velocity.y() - body1.position.y()*body1.velocity.x());

         m_kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();
    }
         m_angularMomentum[0] += L_x;  m_angularMomentum[1] += L_y;  m_angularMomentum[2] += L_z; 
         
}

void SolarSystem::calculateGravitationalForce()
{
    for(CelestialBody &body : m_bodies) {
        // Reset forces on all bodies
        body.force.zeros();
    }

    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body1 = m_bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++) {
            CelestialBody &body2 = m_bodies[j];
            vec3 deltaRVector = body1.position - body2.position;
            double dr = deltaRVector.length();
            vec3 e_r = deltaRVector.normalized();
            // Calculate the force
            body1.force += -G*body1.mass*body2.mass*e_r/pow(dr,2); //beta is an exponent here
            body2.force += G*body1.mass*body2.mass*e_r/pow(dr,2); //beta is an exponent here
        }
       }
         
         
}

void SolarSystem::calculateGravitationalForce_relativistic()
{
    for(CelestialBody &body : m_bodies) {
        // Reset forces on all bodies
        body.force.zeros();
    }
    
    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body1 = m_bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++) {
            CelestialBody &body2 = m_bodies[j];
            vec3 deltaRVector = body1.position - body2.position;
            double dr = deltaRVector.length();
            vec3 e_r = deltaRVector.normalized();
            // Calculate the force
            L_x_Merc = body2.position.y()*body2.velocity.z() - body2.position.z()*body2.velocity.y();
            L_y_Merc = body2.position.z()*body2.velocity.x() - body2.position.x()*body2.velocity.z();
            L_z_Merc = body2.position.x()*body2.velocity.y() - body2.position.y()*body2.velocity.x();
            L_Merc = L_x_Merc*L_x_Merc + L_y_Merc*L_y_Merc + L_z_Merc*L_z_Merc; //l^2
            body2.force += (G*body1.mass*body2.mass*e_r/pow(dr,2))*(1. + (3.*L_Merc)/(dr*dr*c*c)); //keep Sun at rest, so only force on body.2 = Mercury
        }
       }
         
         
}

int SolarSystem::numberOfBodies() const
{
    return m_bodies.size();
}

double SolarSystem::totalEnergy() const
{
    return m_kineticEnergy + m_potentialEnergy;
}

double SolarSystem::potentialEnergy() const
{
    return m_potentialEnergy;
}

double SolarSystem::kineticEnergy() const
{
    return m_kineticEnergy;
}

void SolarSystem::writeToFile(string filename)
{
    if(!m_file.good()) {
        m_file.open(filename.c_str(), ofstream::out);
        if(!m_file.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }

    int body_number = 1;
    for(CelestialBody &body : m_bodies) {
        m_file << body_number << "  " << body.position.x() << " " << body.position.y() << " " << body.position.z() << " " << body.position.length()<< "\n";
     body_number = body_number + 1;
    }
}

vec3 SolarSystem::angularMomentum() const
{
    return m_angularMomentum;
}

std::vector<CelestialBody> &SolarSystem::bodies()
{
    return m_bodies;
}


void SolarSystem::writeToFile_Energy_Momentum(string filename, double initial_total_Energy, double initial_angular_momentum, double i_kin, double i_pot)
{
    if(!n_file.good()) {
        n_file.open(filename.c_str(), ofstream::out);
        if(!n_file.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }   
        
        
        vec3 L = angularMomentum(); //vec3 L = m_angularMomentum;
        double total_Energy = m_kineticEnergy + m_potentialEnergy;  //totalEnergy();
        double L_length;
        L_length = L.length();
        n_file << abs((totalEnergy()-initial_total_Energy)/initial_total_Energy) << " " << abs((L_length-initial_angular_momentum)/initial_angular_momentum) << " " <<m_kineticEnergy<< " " << m_potentialEnergy<<"\n";
 
    
}


class Euler
{
public:
    double m_dt;
    Euler(double dt);
    void integrateOneStep(class SolarSystem &system);
};

Euler::Euler(double dt) :
    m_dt(dt)
{

}

void Euler::integrateOneStep(SolarSystem &system)
{
    system.calculateGravitationalForce();
    int i = 0;
    for(CelestialBody &body : system.bodies()) {
        if (i != 0 ){
        body.position += body.velocity*m_dt;
        body.velocity += body.force / body.mass * m_dt;
       }
       i += 1;
    }
}

class Verlet
{
public:
    double m_dt_Verlet;
    Verlet(double dt);
    void integrateOneStep_Verlet(class SolarSystem &system);
};

Verlet::Verlet(double dt) :
    m_dt_Verlet(dt)
{

}

void Verlet::integrateOneStep_Verlet(SolarSystem &system)
{
    system.calculateGravitationalForce(); //or system.calculateGravitationalForce_relativistic();
    int i = 0; //keep the Sun fixed by omitting the first body
    for(CelestialBody &body : system.bodies()) {
       // if (i != 0 ){
        vec3 a_i = body.force/body.mass;
        body.position +=  body.velocity*m_dt_Verlet + 0.5*m_dt_Verlet*m_dt_Verlet*a_i;
        system.calculateGravitationalForce(); //calculate new acceleration, or system.calculateGravitationalForce_relativistic();
        vec3 a_i_plus_1 = body.force/body.mass;
        body.velocity += 0.5*m_dt_Verlet*(a_i + a_i_plus_1 );
  // }
       i += 1;
}
}


int main()
{
    double years = 150.;
    int numTimesteps = 1.5e6;
    double dt = years/numTimesteps;
   
    SolarSystem solarSystem;
    //CelestialBody &sun = solarSystem.createCelestialBody( vec3(0,0,0), vec3(0,0,0), 1.0 ); //Sun, use for d), e) and g)
    //solarSystem.createCelestialBody( vec3(1, 0, 0), vec3(0,2*M_PI, 0), 3e-6 ); //v_circ = 2*M_PI //Earth, use for d)
        //solarSystem.createCelestialBody( vec3(1.779e-1+3.631e-3, 9.752e-1-7.478e-3, -2.539e-5-1.841e-5), vec3(-1.719e-2, 3.105e-3, -4.95e-7)*365.25, 3e-6 ); //Earth, e)
    //solarSystem.createCelestialBody( vec3(3.738e-1+3.631e-3 , -5.214-7.478e-3 , 1.326e-2-1.841e-5), vec3(7.433e-3,8.997e-4, -1.7005e-4)*365.25, 9.5e-4 ); //Jupiter, mass=9.5e-4, use for e)

    solarSystem.createCelestialBody( vec3(-3.631e-3, 7.478e-3, 1.841e-5), vec3(-8.4e-6, -1.78e-6, 2.315e-7)*365.25, 1.0 ); //Sun moving, NASA -use for f)
    solarSystem.createCelestialBody( vec3(-3.891e-1, -1.63e-1, 2.145e-2), vec3(5.558e-3,-2.451e-2,-2.513e-3)*365.25, 1.65e-7 ); //Mercury
    solarSystem.createCelestialBody( vec3( 6.403e-1,-3.292e-1,-4.1763e-2), vec3(9.245e-3, 1.784e-2,-2.89e-4)*365.25, 2.45e-6 ); //Venus
    solarSystem.createCelestialBody( vec3(1.779e-1, 9.752e-1, -2.539e-5), vec3(-1.719e-2, 3.105e-3, -4.95e-7)*365.25, 3e-6 ); //Earth
    solarSystem.createCelestialBody( vec3(-1.47, -6.58e-1, 2.206e-2), vec3(6.2965e-3, -1.155e-2,-3.964e-4)*365.25, 3.3e-7 ); //Mars
    solarSystem.createCelestialBody( vec3(3.738e-1, -5.214, 1.326e-2), vec3(7.433e-3,8.997e-4, -1.7005e-4)*365.25, 9.5e-4 ); //Jupiter, mass=9.5e-4 
    solarSystem.createCelestialBody( vec3(3.7,-9.322,1.4951e-2), vec3(4.878e-3, 2.04e-3, -2.3e-4)*365.25, 2.75e-4 ); //Saturn
    solarSystem.createCelestialBody( vec3(1.6267e+1, 1.13257e+1,-1.687e-1), vec3(-2.2762e-3, 3.04455e-3, 4.0844e-5)*365.25, 4.4e-5 ); //Uranus
    solarSystem.createCelestialBody( vec3(2.9226e+1,-6.42146,-5.41305e-1), vec3(6.53205e-4, 3.085e-3,-7.855e-5)*365.25, 5.15e-5 ); //Neptune
    solarSystem.createCelestialBody( vec3(1.2913e+1, -3.1373e+1, -3.782e-1), vec3(2.974e-3, 5.2698e-4, -9.16656e-4)*365.25, 6.55e-9 ); //Pluto

    //solarSystem.createCelestialBody( vec3(0.3075, 0, 0), vec3(0,12.44,0), 1.65e-7 ); //Mercury for g)
    vector<CelestialBody> &bodies = solarSystem.bodies();

    for(int i = 0; i<bodies.size(); i++) {
        CelestialBody &body = bodies[i]; // Reference to this body
        cout << "The position of this object is " << body.position << " with velocity " << body.velocity << endl;
    }
    // remember intial values to compare their evolution
    solarSystem.calculateEnergyMomentum();
    double initial_total_Energy = solarSystem.totalEnergy(); 
    vec3 initial_angular_momentum_vec = solarSystem.angularMomentum();
    double initial_angular_momentum = initial_angular_momentum_vec.length();
    double i_kin = solarSystem.potentialEnergy();
    double i_pot = solarSystem.kineticEnergy();
    cout << "Initial angular momentum: " << initial_angular_momentum << endl;
    
    //Here choose algorithm:
    //Euler integrator(dt);
    Verlet integrator_Verlet(dt);
    clock_t start, finish;
    start = clock();
    int counter = 0;
    for(int timestep=0; timestep<numTimesteps; timestep++) {
        integrator_Verlet.integrateOneStep_Verlet(solarSystem); //integrator.integrateOneStep(solarSystem);   
        counter += 1;
        solarSystem.calculateEnergyMomentum();
       // if (counter > 1e8-1e6){ //record only last year for Mercury precession
        solarSystem.writeToFile("positions_model.xyz");
        
//}
        //solarSystem.writeToFile_Energy_Momentum("energy_momentum_Euler_report.txt", initial_total_Energy, initial_angular_momentum, i_kin, i_pot);
    }
    finish = clock();
    double time;
    time = (double(finish - start)/CLOCKS_PER_SEC);
    cout << "Time elapsed: " << time << " s" << endl;

    cout << "The system has " << solarSystem.bodies().size() << " objects." << endl;


    return 0;
}
