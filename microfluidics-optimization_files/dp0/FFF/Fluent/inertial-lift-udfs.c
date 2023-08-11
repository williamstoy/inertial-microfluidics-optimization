/* UDF for computing the lateral displacement of particle suspended in fluid flow */


/**
 * References
 * Poiseuille Flow in Rectangular Channel Calculation (Umax)
 * 	https://web.archive.org/web/20220129140030/http://skinny.dk/picasso_fish/reports/Pedersen%20%282003%29%20Poiseuille%20Flow.pdf
 * Coefficient of Lift and Lift Force Calculations:
 * 	DOI: 10.1039/d1lc00225b
*/

#include "udf.h"
#include "mem.h" // cell indexing header
#include "dpm.h" // particle properties even though declared in the macro arguments
#include "para.h" // parallel calculations
#include <math.h>
#include <stdlib.h>

#include <time.h>

#define BETA(n,H) (2*n-1)*M_PI/H
#define SIGN(x) (int)(x > 0) - (int)(x < 0)
#define REMOVE_PARTICLES FALSE

/** TODO
 * How to pass in the powers and coefficients w/o hard coding
 * Share calculation of Pressure Gradient between UDF functions / files
 * How do we get the domain of the interior of the model?
 * Fix index in lift coefficient calculation. (Will likely need to re-port)
 */

/** USER DEFINED MEMORY LOCATIONS
* 0: CLy
* 1: CLz
* 2: Max X velocity (m/s)
*/


double get_random()
{
  return (double)rand() / RAND_MAX;
}

 
// from http://c-faq.com/lib/gaussian.html, Knuth solution, #3
double gaussrand(double mu, double sigma)
{
  static double V1, V2, S;
  static int phase = 0;
  double X;

  if(phase == 0) {
    do {
      double U1 = get_random();
      double U2 = get_random();

      V1 = 2 * U1 - 1;
      V2 = 2 * U2 - 1;
      S = V1 * V1 + V2 * V2;
      } while(S >= 1 || S == 0);

    X = V1 * sqrt(-2 * log(S) / S);
  } else
    X = V2 * sqrt(-2 * log(S) / S);

  phase = 1 - phase;

  return X * sigma + mu;
}



void CALCULATE_LIFT_COEFFICIENTS(double AR, double Re, double kappa, double p[3], double *CL, int index)
{
  int i, k;
  double neuron_input_sum, sign;
  double inputs[5], outputs_layer_1[20], outputs_layer_2[8];

  static const double offsets_layer_1[20] = { -3.0912897985545538, 3.3260267431311772,
    2.5103508876609535, -1.6013936251390251, -1.6577405176451285,
    -1.4676602992259893, 0.79308293598130486, 0.38060470412883013,
    0.052627557434173376, -0.11942570835094345, -0.10627757338094153,
    0.21447627159805541, 0.55768068330094722, -1.300580094035618,
    -0.60919426751306238, -1.773497218413272, 0.32350624950536855,
    0.031700372409363284, -4.27049712836732, 1.5309232034203801 };

  static const double weights_layer_1[100] = { 0.4629377063578074, -0.41713620529734652,
    0.033754762283431591, 0.41339623376264856, 0.69791559374472112,
    0.079269760502799172, 0.23130377988890027, 0.0068065268541734978,
    1.5360717676570379, -1.1493762140461834, 0.031803192942887422,
    -0.14274969007125285, 0.21185146045540321, -0.33162851183309666,
    -0.0427928083054695, -1.0617927846762474, 1.0417110746023555,
    -0.16679670035121544, -5.2899410860527016, 0.074018728390475755,
    -0.096596165030048461, -0.07061134172387247, -0.18467142405730075,
    -0.035631726727613805, 0.073732880418435287, -0.027325176815033141,
    0.14493179276776816, 0.034553309012998842, -0.097496480492810192,
    -0.10647994128154212, -0.19605480358714023, 0.00082205616236863594,
    0.0618845229700362, -0.038624577837312678, -0.26499828121862257,
    0.57212512997881571, -0.092768766623457483, 0.062875456906280158,
    -0.49444476278653132, -0.29646281033588884, -0.0659914675489781,
    0.022783355786681534, -0.028425522112777636, -0.095243945772698846,
    -0.056894965632009645, 0.15717407278255258, 0.071901533343011442,
    -0.0084678146126629175, -0.069745921198961988, 0.031688610296034515,
    0.23893011428979116, 0.0067908671985627729, -0.016432387720108917,
    -0.0047202951710107873, -0.17999436073129729, 0.18053647734515868,
    -0.024547913645244872, -0.10169838832446142, 0.19092927332882839,
    0.031512531652647634, -0.24377336846899061, -3.132767759887193,
    -0.12772629452475445, -0.797267944895747, 1.790627772819086,
    -0.21820130204417845, 0.22937576408324722, -0.3671065440154927,
    -0.3235991332391141, -1.2191295750274629, 0.80177870126230322,
    -0.2224638682511173, 0.28140456528490954, -0.85971750101531552,
    0.37202547008967285, 0.055121375470426628, -0.36590477070650923,
    1.8330694288866971, 0.28467917229722506, -0.76695829812169,
    2.0163005537310448, 0.17186555038230611, -1.7495943026222942,
    -0.90615304994873569, -1.1239634037792847, -1.5834329186740805,
    0.7641015918262557, -0.51093222468839139, 0.038604050024509638,
    0.58203063030696722, 0.058960801322799745, 0.54957139932298027,
    -0.24583412187077791, -0.37254646382630924, 0.28961127364562883,
    -0.33800776725353521, -0.17462268555262878, 0.25170262247740482,
    -0.60521930361883947, -0.069477830002075441 };

  static const double offsets_layer_2[8] = { -0.98247358224342662, 1.6059512044901776,
    -0.009036195639525315, -0.37216527387034126, -2.3442114361903883,
    -3.273819351472957, 2.3833600596693563, 0.46610841161464789 };

  static const double weights_layer_2[160] = { 1.1747672816166976, 0.60336888706423386,
    1.2867669361124812, -0.34468711708390259, -1.3968612722553262,
    -0.64000077114090037, -0.73132466211602332, -0.053305699218618538,
    -0.48205135999910487, 0.1055348913703812, -1.1339464708998355,
    0.19291689507467269, -0.28684366576348874, -0.27365377350330794,
    -0.1839120291402663, 0.21580034629065659, -2.3799689055011992,
    1.0853416102936644, -1.3899584455425995, -0.99977823083062023,
    0.7840552088438455, 1.4239399412999303, -1.564772967056383,
    0.86622318097537609, -1.903204521819339, 1.0360697523701035,
    -1.1240286053829236, 1.0990817930486285, 0.53642939668499234,
    0.14257478957425465, -0.60524680013375187, 0.23744657696078608,
    -1.088076896412757, 0.64145105656798029, -2.2523027577080676,
    -0.18307838624024572, 0.07280446293692977, -0.099802273329160821,
    -0.75825413776726314, 0.25989808268870973, -1.7248482528290132,
    0.48978362907026296, -0.40921528724012113, -0.57997352600178431,
    -0.1787467255444167, -0.16341667367539742, -0.5527483319594626,
    0.1152982170189759, -1.2995350997416959, 0.18609684769264048,
    -0.45382010310634696, -0.72497198362425519, -0.34389428489547774,
    -0.5013665050229581, -0.15426602173114456, 1.6781945490731887,
    -2.1528859421055762, 1.3495421152219236, -1.7198269606235717,
    -1.1181947103407623, 0.94776506564412522, 0.57764322960493319,
    -3.8666098074407191, 1.5833372851588574, -0.29605969015368339,
    -0.376203580073117, -0.73869560466901985, 0.16663062691129243,
    -0.32219740402085911, -0.20646892346321052, -0.21223344512062853,
    0.75760940432625012, -0.0068295042610621633, -0.11972160121134233,
    0.14600085710456995, 0.084642727314444069, -0.25819709122096074,
    0.0879476723195499, 0.056503558600723786, 0.017598701380658752,
    -0.16228296619835317, -0.30726238607751483, -0.22617408954115162,
    0.33286475050470565, 1.4334731218177841, 0.35522684846591707,
    -0.31887437741816449, 0.20136696633496534, -0.16241008517167643,
    0.39641922804642943, 0.045006513795005622, 0.64230055261714192,
    0.658062186632713, 0.64082615370307006, -4.0828682474327325,
    -1.6215902898916303, -3.1766921348543997, -0.38066138769783342,
    -0.026397366538768768, 0.573718639240188, 0.13937169624292156,
    1.6802973480227799, -3.5367811476006654, 1.0131448722834098,
    -0.98826075804864688, -1.3042893409897465, -0.47761827745626273,
    -0.17560872016815784, 0.12956003513819755, -0.36712028434308208,
    0.72220354455444058, 1.1841540940352397, -0.15653848012037783,
    -0.51466323522196855, -0.070617884440168185, 0.56176219770192415,
    0.94234708921382759, 0.938890346670385, 0.180232710113337,
    -0.0067430807963181965, 0.68285627607651478, 0.010387364350044054,
    0.1137418523711472, -0.0135938438311762, 0.045982827980147981,
    0.057169961999039579, 0.35917462047381016, -0.48129418730721413,
    0.48974693540572395, 0.73252226354749894, 1.180505997849048,
    -0.48491657154843854, 0.51307137073346709, 0.3948887129206719,
    0.25881084558793471, -1.2192428623326246, -0.050254992618655334,
    0.048873559346412306, 0.063179651198612491, -0.14847429881066151,
    0.28826375900604995, -0.28080203741735926, -0.13994694414804182,
    0.014561437916550395, 0.34157781125919734, -0.001356881569425151,
    -0.1489353919577526, 0.040500905971888851, 0.11910591096808623,
    0.13192138581886556, 0.14324657587058418, -0.28855608132209931,
    0.94271449144581576, -0.45004710846411133, 0.061179451165720246,
    0.28427875312043527, -1.4464009437863961, -0.84024239337213713,
    0.4370360717375803, -0.14302306101141724 };

  static const double weights_layer_3[16] = { 0.18826376956017257, -3.3997287111751024,
    -2.5194013888960471, -1.9826222806895024, 0.19009832278140057,
    -0.17292280720160841, -2.3843096996956308, -0.063198924695663933,
    -4.5685103708787951, 0.026365975274439686, 4.0760256849122589,
    0.016385163105698302, 0.054440678326434823, -1.7795902634438079,
    0.570248458711622, -4.3548445693126441 };

  /*  ===== INPUT NORMALIZATION ===== */
  inputs[0] = (AR - 1.0) * 0.33333333333333331 * 2.0 - 1.0;
  inputs[1] = (Re - 50.0) * 0.0066666666666666671 * 2.0 - 1.0;
  inputs[2] = (kappa - 0.1) * 5.0 * 2.0 - 1.0;
  inputs[3] = fabs(p[1]) * 1.25 * 2.0 - 1.0;
  inputs[4] = fabs(p[2]) * 1.0526315789473679 * 2.0 - 1.0;

  /*  Layer 1 */
  for (k = 0; k < 20; k++) {
    neuron_input_sum = 0.0;
    for (i = 0; i < 5; i++) {
      neuron_input_sum += weights_layer_1[k + 20 * i] * inputs[i];
    }

    outputs_layer_1[k] = tanh(offsets_layer_1[k] + neuron_input_sum);
  }

  /*  Layer 2 */
  for (k = 0; k < 8; k++) {
    neuron_input_sum = 0.0;
    for (i = 0; i < 20; i++) {
      neuron_input_sum += weights_layer_2[k + (i << 3)] * outputs_layer_1[i];
    }

    outputs_layer_2[k] = tanh(offsets_layer_2[k] + neuron_input_sum);
  }

  /*  Layer 3 */
  neuron_input_sum = 0.0;
  for (k = 0; k < 8; k++) {
    neuron_input_sum += weights_layer_3[(index-1) + (k << 1)] * outputs_layer_2[k];
  }

  sign = 1;
  if (p[index] < 0) {
    sign = -1;
  }

  /*  ===== OUTPUT DENORMALIZATION ===== */
  *CL = sign * ((((1.6035102046923548 * (double)(index-1) + -0.13413896296737876) + neuron_input_sum) + 1.0)
    / 2.0 / (0.093013622399621965 * (double)(index-1) + 0.99451030312674) +
    (0.085999999999999965 * (double)(index-1) + -0.7776));
}



double CALCULATE_PRESSURE_GRADIENT(double W, double H, double mu, double Q)
{
  // DECLARATION OF VARIABLES
  double dPdx, s, k;
  int n;

  // calculate the pressure gradient at this point given the flow rate
  s = 0;
  for (n = 1; n <= 5; n++) {
      s += (1/pow(2*n-1,5))* (cosh(BETA(n,H) * W) - 1) / sinh(BETA(n,H) * W);
  }
  k = ((pow(H,3))*W)/(12*mu) - (16*pow(H,4))/(pow(M_PI,5)*mu) * s;
  dPdx = Q/k; // Pressure gradient

  return dPdx;
}


/***
 *  Solution to PDEs of 3D Poiseuille Flow in a Rectangular Channel
 *  Boussinesq derivation of flow velocity in rectangular channel
 *  
 * 	arguments:
 * 	double w: channel width in meters
 *  double h: channel height in meters
 *  double Q: volumetric flow rate in the channel
 */
double CALCULATE_MAX_VELOCITY(double H, double W, double mu, double dPdx)
{
	// DECLARATION OF VARIABLES
	double Umax, s;
	int n;

	s = 0;
	for (n = 1; n <= 5; n++) {
	    s += (1/pow(2*n-1,3)) * ((2*sinh(BETA(n,H) * W / 2)) / sinh(BETA(n,H) * W));
	}

	Umax = (dPdx/(2*mu))*(pow(H/2,2)) - ((4*dPdx*pow(H,2))/(mu*pow(M_PI,3))) * s;

	// return the velocity in the x direction
	return Umax;
}


double CALCULATE_VELOCITY(double H, double W, double y, double z, double mu, double dPdx)
{
	// DECLARATION OF VARIABLES
	double unm, ux;
	int n, m, maxiter;

	unm = 0;
	maxiter = 20;
	maxiter = (maxiter + (1 - maxiter % 2));
	for (n = 1; n <= maxiter; n += 2) {
		for (m = 1; m <= maxiter; m += 2) {
			unm += (1 / (n * m * (pow(n, 2) / pow(W, 2) + pow(m, 2) / pow(H, 2)))) * sin(n * (M_PI / W) * z) * sin(m * (M_PI / H) * y);
		}
	}
	ux = (16 / pow(M_PI, 4)) * (dPdx / mu) * unm;

	return ux;
}


void CALCULATE_CHANNEL_PARAMETERS(double particle_diameter, double mu, double rho, double *H, double *W, double *Re, double *kappa, double particle_position[3], double channel_center_offset[3])
{
	double x = particle_position[0];
	double A, P, Q, Dh, yoffset, zoffset;

	Q = Get_Input_Parameter("volumetric_flow_rate_ul_per_min"); // uL/min
  //Message("Volumetric Flow Rate (User Parameter): %.2f [uL/min]\n", Q);
  Q = Q * 1e-9 / 60; // convert Q to m^3/s

	/* THIS IS HOW YOU DESCRIBE THE CHANNEL AR */
	double notch_length = Get_Input_Parameter("notch_length");
  double notch_height = Get_Input_Parameter("notch_height");
	double notch_midpoint = Get_Input_Parameter("notch_midpoint");
  double channel_width = Get_Input_Parameter("channel_width");
  double channel_height = Get_Input_Parameter("channel_height");

  *W = channel_width;
  *H = channel_height;
  yoffset = 0;
  zoffset = 0;

	if (x > (notch_midpoint - notch_length/2) && x <= (notch_midpoint + notch_length/2)) { // we are in the notch
    // we need to shift the center offset down due to the presence of the notch
    yoffset = notch_height/2;
    zoffset = 0;

		*H = *H - notch_height;
	}

  channel_center_offset[0] = 0.0;
  channel_center_offset[1] = yoffset;
  channel_center_offset[2] = zoffset;

	/* normalize the particle position in y and z with respect to the center position of the channel */
	particle_position[1] = 2 * (particle_position[1] + yoffset) / *H;
	particle_position[2] = 2 * (particle_position[2] + zoffset) / *W;

	A = *W * *H; // channel area
	P = *W * 2 + *H * 2; // channel perimeter

	Dh = 4 * A / P; // hydraulic
	*Re = rho * Q * Dh / (mu * A);
	*kappa = particle_diameter / *H; // blockage ratio
}


DEFINE_DPM_BODY_FORCE(inertial_lift, p, i)
{
  if (i == 0) {
    return 0.0; // do not calculate lift force in the X direction
  }

	// DECLARATION OF VARIABLES
	double CL, FL, particle_diameter, H, W, Re, kappa, mu, rho, Umax, x, y, z;
	cell_t c = P_CELL(p); // the cell in which the particle is present
	Thread *t = P_CELL_THREAD(p); // thread initialization

  double channel_center_offset[3];
  channel_center_offset[0] = 0;
  channel_center_offset[1] = 0;
  channel_center_offset[2] = 0;

	double particle_position[3];
  x = P_POS(p)[0]; y = P_POS(p)[1]; z = P_POS(p)[2];
  particle_position[0] = x;
  particle_position[1] = y;
  particle_position[2] = z;

	particle_diameter = P_DIAM(p); 		// Particle diameter
	mu = C_MU_L(c,t); 	// Cell dynamic viscosity
	rho = C_R(c,t); 	// Cell density

	CALCULATE_CHANNEL_PARAMETERS(particle_diameter, mu, rho, &H, &W, &Re, &kappa, particle_position, channel_center_offset);

  // figure out if we should remove this particle from the list (it has hit the wall)
  double notch_midpoint, notch_length;
  notch_midpoint = Get_Input_Parameter("notch_midpoint");
  notch_length = Get_Input_Parameter("notch_length");

  double channel_height = Get_Input_Parameter("channel_height");
  double channel_width = Get_Input_Parameter("channel_width");
  double notch_height = Get_Input_Parameter("notch_height");

  double zmax = (channel_width/2 - particle_diameter/2 - channel_center_offset[2]);
  double zmin = (-channel_width/2 + particle_diameter/2 - channel_center_offset[2]);
  double ymax = (channel_height/2 - particle_diameter/2 - channel_center_offset[1]);
  double ymin = (-channel_height/2 + particle_diameter/2 - channel_center_offset[1]);

  if (y < ymin) {
    // we have hit a top or bottom wall, reflect the y velocity
    P_VEL(p)[1] = abs(P_VEL(p)[1]);
    P_POS(p)[1] = ymin + 1e-7; // reset the position inside the channel but inset by 100nm
  } else if (y > ymax) {
    // we have hit a top or bottom wall, reflect the y velocity
    P_VEL(p)[1] = -abs(P_VEL(p)[1]);
    P_POS(p)[1] = ymax - 1e-7; // reset the position inside the channel but inset by 100nm
  }

  if (z < zmin) {
    // we have hit a side wall, reflect the z velocity
    P_VEL(p)[2] = abs(P_VEL(p)[2]);
    P_POS(p)[2] = zmin + 1e-7; // reset the position inside the channel but inset by 100 nm
  } else if (z > zmax) {
    // we have hit a side wall, reflect the z velocity
    P_VEL(p)[2] = -abs(P_VEL(p)[2]);
    P_POS(p)[2] = zmax - 1e-7; // reset the position inside the channel but inset by 100 nm
  }

  double notch_leading_edge_x = notch_midpoint - notch_length/2;
  double notch_bottom_edge_y = channel_height/2 - notch_height;
  // calculate what happens at the corner
  if (x >= notch_leading_edge_x - particle_diameter/2 && x < notch_leading_edge_x && y >= notch_bottom_edge_y) {
    P_VEL(p)[0] = -abs(P_VEL(p)[0]);
    P_POS(p)[0] = notch_leading_edge_x - particle_diameter/2 - 1e-7;
  }

  double distance_to_leading_edge = sqrt(pow(P_POS(p)[0] - notch_leading_edge_x, 2) + pow(P_POS(p)[1] - notch_bottom_edge_y, 2));

  if(distance_to_leading_edge <= particle_diameter/2 && P_POS(p)[0] < notch_leading_edge_x && P_POS(p)[0] > notch_bottom_edge_y) {
    double q, vx, cx, rx, vy, cy, ry;
    vx = P_VEL(p)[0];
    vy = P_VEL(p)[1];
    cx = P_POS(p)[0];
    cy = P_POS(p)[1];
    rx = notch_leading_edge_x;
    ry = notch_bottom_edge_y;
    q = -2 * (vx*(cx-rx) + vy*(cy-ry)) / pow(particle_diameter/2, 2);
    P_VEL(p)[0] = vx + q * (cx - rx);
    P_VEL(p)[1] = vy + q * (cy - ry);
  }

	CALCULATE_LIFT_COEFFICIENTS(W/H, Re, kappa, particle_position, &CL, i);
  /* Potential time savings using precomputed lift coefficient values stored in cells
   * In a quick test, this did not save a significant amount of time
   */
  //CL = C_UDMI(c, t, i-1);

	// Maximum velocity for each cell should be precalculated prior to starting DPM calculation
	Umax = C_UDMI(c, t, 2);

	// Calculate the lift force based on Su et al 2021 
	FL = CL * rho * pow(Umax, 2) * pow(particle_diameter, 4) / pow(H, 2);

	// An acceleration should be returned
	return (FL / P_MASS(p));
}


DEFINE_PROFILE(inlet_x_velocity, thread, position)
{
	// DECLARATION OF VARIABLES
	double face_centroid_pos[3]; // this will hold the position vector
	double dPdx, W, H, y, z, mu, Q;
	face_t f;

	W = Get_Input_Parameter("channel_width");
	H = Get_Input_Parameter("channel_height");

	// TODO: would be nice to not have to hard code this. Can get this from the cell, but 
	mu = 0.001003; // kg/(m*s)
	Q = Get_Input_Parameter("volumetric_flow_rate_ul_per_min"); // uL/min
  //Message("Volumetric Flow Rate (User Parameter): %.2f [uL/min]\n", Q);
  Q = Q * 1e-9 / 60; // convert Q to m^3/s

	dPdx = CALCULATE_PRESSURE_GRADIENT(H, W, mu, Q);

	begin_f_loop(f, thread)
	{
		F_CENTROID(face_centroid_pos, f, thread);
		
		y = face_centroid_pos[1] + 0.5 * H;
		z = face_centroid_pos[2] + 0.5 * W;
		F_PROFILE(f, thread, position) = CALCULATE_VELOCITY(H, W, y, z, mu, dPdx);
	}
	end_f_loop(f, thread)

	Message("Finished Loading Profile\n");
}



DEFINE_INIT(precalculate_inertial_lift_coefficients_and_umax, domain)
{
	// DECLARATION OF VARIABLES
	double cell_centroid_position[ND_ND]; /* this will hold the position vector */
	double dPdx, W, H, x, y, z, mu, Q, Umax, particle_diameter, rho, Re, kappa, CL;
	face_t f;
	cell_t c;
	Thread *t;

  double channel_center_offset[3];
  channel_center_offset[0] = 0;
  channel_center_offset[1] = 0;
  channel_center_offset[2] = 0;

	Q = Get_Input_Parameter("volumetric_flow_rate_ul_per_min"); // uL/min
	Message("Volumetric Flow Rate (User Parameter): %.2f [uL/min]\n", Q);
	Q = Q * 1e-9 / 60; // convert Q to m^3/s

  particle_diameter = Get_Input_Parameter("particle_diameter");

  Message("Particle Diameter (User Parameter): %.2f [um]\n", particle_diameter * 1e6);

	thread_loop_c(t, domain)
	{
		begin_c_loop(c, t)
		{
			C_CENTROID(cell_centroid_position, c, t); // get the cell position
			mu = C_MU_L(c,t); 	// Cell dynamic viscosity
			rho = C_R(c,t); 	// Cell density

			CALCULATE_CHANNEL_PARAMETERS(particle_diameter, mu, rho, &H, &W, &Re, &kappa, cell_centroid_position, channel_center_offset);

      // from the height and width of the channel, determine what the pressure gradient and maximum velocity are
      dPdx = CALCULATE_PRESSURE_GRADIENT(H, W, mu, Q);
      Umax = CALCULATE_MAX_VELOCITY(H, W, mu, dPdx);

      // calculate the coefficient of lift for a particle at the cell center
      CALCULATE_LIFT_COEFFICIENTS(W/H, Re, kappa, cell_centroid_position, &CL, 1);
      C_UDMI(c, t, 0) = CL;
      CALCULATE_LIFT_COEFFICIENTS(W/H, Re, kappa, cell_centroid_position, &CL, 2);
      C_UDMI(c, t, 1) = CL;

			// save the maximum velocity of the channel at this location
			C_UDMI(c, t, 2) = Umax;
		}
		end_c_loop(c, t)
	}

	Message("Finished Setting Umax for all Cells in Domain\n");

  // seed the RNG (important for later)
  srand(time(NULL));
}



DEFINE_ON_DEMAND(port_input_vars_to_scheme)
{
  // port the input parameters to scheme parameters
  double channel_width, channel_height, channel_length, notch_height, notch_length, re, max_number_of_fluent_iterations, number_of_steps, volumetric_flow_rate_ul_per_min;
  
  channel_width = Get_Input_Parameter("channel_width");
  channel_height = Get_Input_Parameter("channel_height");
  channel_length = Get_Input_Parameter("channel_length");
  notch_height = Get_Input_Parameter("notch_height");
  notch_length = Get_Input_Parameter("notch_length");
  re = Get_Input_Parameter("re");
  max_number_of_fluent_iterations = Get_Input_Parameter("max_number_of_fluent_iterations");
  number_of_steps = Get_Input_Parameter("number_of_steps");
  volumetric_flow_rate_ul_per_min = Get_Input_Parameter("volumetric_flow_rate_ul_per_min");

  RP_Set_Real("channel_width_scheme", (real) (1e6 * channel_width));
  RP_Set_Real("channel_height_scheme", (real) (1e6 * channel_height));
  RP_Set_Real("channel_length_scheme", (real) (1e6 * channel_length));
  RP_Set_Real("notch_height_scheme", (real) (1e6 * notch_height));
  RP_Set_Real("notch_length_scheme", (real) (1e6 * notch_length));
  RP_Set_Real("re_scheme", (real) re);
  RP_Set_Real("max_number_of_fluent_iterations_scheme", (real) max_number_of_fluent_iterations);
  RP_Set_Real("number_of_steps_scheme", (real) number_of_steps);
  RP_Set_Real("volumetric_flow_rate_ul_per_min_scheme", (real) volumetric_flow_rate_ul_per_min);

  Message("Ported Input Variables to Scheme\n");
}



DEFINE_DPM_INJECTION_INIT(particle_init,I)
{
  double H, W, mu, Q, dPdx, d, d2, y, z, r, s;
  Particle *p;
  int count = 0;
  int max_count = Get_Input_Parameter("max_dpm_particles");

  srand(time(NULL)); // randomize seed

  W = Get_Input_Parameter("channel_width");
  H = Get_Input_Parameter("channel_height");

  mu = 0.001003; // kg/(m*s)
  Q = Get_Input_Parameter("volumetric_flow_rate_ul_per_min"); // uL/min
  //Message("Volumetric Flow Rate (User Parameter): %.2f [uL/min]\n", Q);
  Q = Q * 1e-9 / 60; // convert Q to m^3/s

  dPdx = CALCULATE_PRESSURE_GRADIENT(H, W, mu, Q);

  double particle_mean = Get_Input_Parameter("particle_diameter");
  double particle_std = Get_Input_Parameter("particle_diameter_std");

  // NB: INSERTING particles on the side 1/3rds of the channel greatly increases the sorting speed
  // and lowers the number of notch features that are required to focus the particles
  // realized this because the particles that focus the slowest start in the top center quadrant
  // set zexclude to 3 to restrict particle insertion to the outer 1/3rds of the channel
  // set zexclude to 1 to put particles across the whole channel inlet width
  double zexclude = 1;

  Message("Initializing Injection: %s\n",I->name);
  loop(p,I->p)  // Standard ANSYS FLUENT Looping Macro to get particle streams in an Injection
  {
    count += 1; // increment the counter

    if (count > max_count) {
      p->stream_index=-1;
      P_VEL(p)[0] = 0;
      P_POS(p)[0] = -100e-6;
    } else {
      // randomize the size of particles
      d = gaussrand(particle_mean, particle_std);
      P_DIAM(p) = d;

      // define particle z locations (across the width of the channel)
      r = get_random();
      y = (r * (H - d) + d/2 ) - H/2;
      r = get_random();
      double r2 = get_random();
      s = (double)SIGN(2 * r2 - 1);

      if (zexclude != 1) {
        z = s * (r * (W/zexclude - d) + (W/(2*zexclude) + d/2));
      } else {
        z = (r * (W - d) + d/2 ) - W/2;
      }
      

      P_POS(p)[1] = y;
      P_POS(p)[2] = z;
      // set the velocity in x based on the particle position
      P_VEL(p)[0] = CALCULATE_VELOCITY(H, W, y+H/2, z+W/2, mu, dPdx);
    }
    
  }
}










