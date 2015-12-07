#include <iostream>
#include "capd/capdlib.h"
using namespace capd;
using namespace std;


// var:lambda_0_0, p_0_0, t_0_0, theta_0_0, thetaI_0_0, w_0_0;par:fc_0_0, i_0_0, pe_0_0, tau_0_0;fun:(4*(((0.9*((((-0.366)+((0.08979*w_0_0)*p_0_0))+((((-0.0337)*w_0_0)*p_0_0)*p_0_0))+(((0.0001*w_0_0)*w_0_0)*p_0_0)))/(1*fc_0_0))+((-1)*lambda_0_0))), (0.41328*(((2*(((2.821+((-0.05231)*theta_0_0))+((0.10299*theta_0_0)*theta_0_0))+((((-0.00063)*theta_0_0)*theta_0_0)*theta_0_0)))*((((p_0_0/1)+((-1)*(((p_0_0/1))^(2)))))^(0.5)))+((-1)*(0.9*((((-0.366)+((0.08979*w_0_0)*p_0_0))+((((-0.0337)*w_0_0)*p_0_0)*p_0_0))+(((0.0001*w_0_0)*w_0_0)*p_0_0)))))), 1, (10*(thetaI_0_0+((-1)*theta_0_0))), (100.000000000000000000000000000000*cos((10*t_0_0))), (50.000000000000000000000000000000*cos((10*t_0_0)));
// set_param: fc_0_0 ==> [0.6537,0.6537]
// set_param: i_0_0 ==> [0,0]
// set_param: pe_0_0 ==> [0,0]
// set_param: tau_0_0 ==> [0,0]
// X_0 : {[14.7,14.7],[0.9833,0.9833],[0,0],[8.8,8.8],[80,80],[110,110]}
// X_t : {[0,100],[0,10],[0,10],[70,70],[0,180],[0,150]}
// T   : [0,10]

int main(){
    cout.precision(17);
    // define an instance of class IMap that describes the vector field
    IMap pendulum("var:lambda_0_0, p_0_0, t_0_0, theta_0_0, thetaI_0_0, w_0_0;par:fc_0_0, i_0_0, pe_0_0, tau_0_0;fun:(4*(((0.9*((((-0.366)+((0.08979*w_0_0)*p_0_0))+((((-0.0337)*w_0_0)*p_0_0)*p_0_0))+(((0.0001*w_0_0)*w_0_0)*p_0_0)))/(1*fc_0_0))+((-1)*lambda_0_0))), (0.41328*(((2*(((2.821+((-0.05231)*theta_0_0))+((0.10299*theta_0_0)*theta_0_0))+((((-0.00063)*theta_0_0)*theta_0_0)*theta_0_0)))*((((p_0_0/1)+((-1)*(((p_0_0/1))^(2)))))^(0.5)))+((-1)*(0.9*((((-0.366)+((0.08979*w_0_0)*p_0_0))+((((-0.0337)*w_0_0)*p_0_0)*p_0_0))+(((0.0001*w_0_0)*w_0_0)*p_0_0)))))), 1, (10*(thetaI_0_0+((-1)*theta_0_0))), (100.000000000000000000000000000000*cos((10*t_0_0))), (50.000000000000000000000000000000*cos((10*t_0_0)));");
    pendulum.setParameter("fc_0_0", 0.6537);
    pendulum.setParameter("i_0_0", 0);
    pendulum.setParameter("pe_0_0", 0);
    pendulum.setParameter("tau_0_0", 0);
    int order=20;
    IOdeSolver solver(pendulum,order);
    ITimeMap timeMap(solver);
    // specify initial condition
    IVector u(6);
    u[0] = 14.7;
    u[1] = 0.9833;
    u[2] = 0.0;
    u[3] = 8.8;
    u[4] = 80;
    u[5] = 110;
    double initTime = 0;
    double finalTime = 1;
    C0Rect2Set set(u,initTime);
    // define functional object
    cerr << "before solution";
    ITimeMap::SolutionCurve solution(initTime);
    cerr << "after solution";
    // and integrate
    timeMap(finalTime,set,solution);
    cerr << "after integrate";
    // then you can ask for any intermediate point on the trajectory
    cout << "domain = [" << solution.getLeftDomain() << "," << solution.getRightDomain() << "]\n";
    cout << "solution(0.5) =" << solution(0.5) << endl;
    // cout << "solution(4) =" << solution(4) << endl;
    // cout << "solution(5.,5.125) = " << solution(interval(5,5.125)) << endl;
    // cout << "solution(7,8) =" << solution(interval(7,8)) << endl;
    // cout << "solution(10)=" << solution(10) << endl;
    return 0;
}
