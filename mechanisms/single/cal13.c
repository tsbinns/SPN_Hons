/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__cal13
#define _nrn_initial _nrn_initial__cal13
#define nrn_cur _nrn_cur__cal13
#define _nrn_current _nrn_current__cal13
#define nrn_jacob _nrn_jacob__cal13
#define nrn_state _nrn_state__cal13
#define _net_receive _net_receive__cal13 
#define rates rates__cal13 
#define states states__cal13 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define pbar _p[0]
#define damod _p[1]
#define maxMod _p[2]
#define level _p[3]
#define max2 _p[4]
#define lev2 _p[5]
#define ical _p[6]
#define m _p[7]
#define h _p[8]
#define ecal _p[9]
#define cali _p[10]
#define calo _p[11]
#define minf _p[12]
#define mtau _p[13]
#define hinf _p[14]
#define htau _p[15]
#define Dm _p[16]
#define Dh _p[17]
#define v _p[18]
#define _g _p[19]
#define _ion_cali	*_ppvar[0]._pval
#define _ion_calo	*_ppvar[1]._pval
#define _ion_ical	*_ppvar[2]._pval
#define _ion_dicaldv	*_ppvar[3]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_ghk(void);
 static void _hoc_modulation(void);
 static void _hoc_rates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_cal13", _hoc_setdata,
 "ghk_cal13", _hoc_ghk,
 "modulation_cal13", _hoc_modulation,
 "rates_cal13", _hoc_rates,
 0, 0
};
#define ghk ghk_cal13
#define modulation modulation_cal13
 extern double ghk( _threadargsprotocomma_ double , double , double );
 extern double modulation( _threadargsproto_ );
 /* declare global and static user variables */
#define q q_cal13
 double q = 2;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "pbar_cal13", "cm/s",
 "ical_cal13", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "q_cal13", &q_cal13,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[4]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"cal13",
 "pbar_cal13",
 "damod_cal13",
 "maxMod_cal13",
 "level_cal13",
 "max2_cal13",
 "lev2_cal13",
 0,
 "ical_cal13",
 0,
 "m_cal13",
 "h_cal13",
 0,
 0};
 static Symbol* _cal_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 20, _prop);
 	/*initialize range parameters*/
 	pbar = 0;
 	damod = 0;
 	maxMod = 1;
 	level = 0;
 	max2 = 1;
 	lev2 = 0;
 	_prop->param = _p;
 	_prop->param_size = 20;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_cal_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cali */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* calo */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ical */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dicaldv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _cal13_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("cal", 2.0);
 	_cal_sym = hoc_lookup("cal_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 20, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 cal13 C:/Users/tomth/SPN_Hons/mechanisms/single/cal13.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96485.3;
 static double R = 8.3145;
static int _reset;
static char *modelname = "LVA L-type calcium current (Cav1.3)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsproto_);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargs_ ) ;
   Dm = ( minf - m ) / mtau * q ;
   Dh = ( hinf - h ) / htau * q ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rates ( _threadargs_ ) ;
 Dm = Dm  / (1. - dt*( ( ( ( ( - 1.0 ) ) ) / mtau )*( q ) )) ;
 Dh = Dh  / (1. - dt*( ( ( ( ( - 1.0 ) ) ) / htau )*( q ) )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rates ( _threadargs_ ) ;
    m = m + (1. - exp(dt*(( ( ( ( - 1.0 ) ) ) / mtau )*( q ))))*(- ( ( ( ( minf ) ) / mtau )*( q ) ) / ( ( ( ( ( - 1.0 ) ) ) / mtau )*( q ) ) - m) ;
    h = h + (1. - exp(dt*(( ( ( ( - 1.0 ) ) ) / htau )*( q ))))*(- ( ( ( ( hinf ) ) / htau )*( q ) ) / ( ( ( ( ( - 1.0 ) ) ) / htau )*( q ) ) - h) ;
   }
  return 0;
}
 
static int  rates ( _threadargsproto_ ) {
    minf = 1.0 / ( 1.0 + exp ( ( v - ( - 33.0 ) ) / ( - 6.7 ) ) ) ;
   mtau = 0.06 + 1.0 / ( exp ( ( v - 10.0 ) / 20.0 ) + exp ( ( v - ( - 17.0 ) ) / - 48.0 ) ) ;
   hinf = 1.0 / ( 1.0 + exp ( ( v - ( - 13.4 ) ) / 11.9 ) ) ;
   htau = 44.3 ;
     return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
double ghk ( _threadargsprotocomma_ double _lv , double _lci , double _lco ) {
   double _lghk;
 double _lz , _leci , _leco ;
 _lz = ( 1e-3 ) * 2.0 * FARADAY * _lv / ( R * ( celsius + 273.15 ) ) ;
   if ( _lz  == 0.0 ) {
     _lz = _lz + 1e-6 ;
     }
   _leco = _lco * ( _lz ) / ( exp ( _lz ) - 1.0 ) ;
   _leci = _lci * ( - _lz ) / ( exp ( - _lz ) - 1.0 ) ;
   _lghk = ( 1e-3 ) * 2.0 * FARADAY * ( _leci - _leco ) ;
   
return _lghk;
 }
 
static void _hoc_ghk(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  ghk ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) );
 hoc_retpushx(_r);
}
 
double modulation ( _threadargsproto_ ) {
   double _lmodulation;
 _lmodulation = 1.0 + damod * ( ( maxMod - 1.0 ) * level + ( max2 - 1.0 ) * lev2 ) ;
   if ( _lmodulation < 0.0 ) {
     _lmodulation = 0.0 ;
     }
   
return _lmodulation;
 }
 
static void _hoc_modulation(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  modulation ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cali = _ion_cali;
  calo = _ion_calo;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cali = _ion_cali;
  calo = _ion_calo;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_cal_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_cal_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_cal_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_cal_sym, _ppvar, 3, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  m = m0;
 {
   rates ( _threadargs_ ) ;
   m = minf ;
   h = hinf ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  cali = _ion_cali;
  calo = _ion_calo;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ical = pbar * m * m * h * ghk ( _threadargscomma_ v , cali , calo ) * modulation ( _threadargs_ ) ;
   }
 _current += ical;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  cali = _ion_cali;
  calo = _ion_calo;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dical;
  _dical = ical;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicaldv += (_dical - ical)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ical += ical ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  cali = _ion_cali;
  calo = _ion_calo;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "cal13.mod";
static const char* nmodl_file_text = 
  "TITLE LVA L-type calcium current (Cav1.3)\n"
  "\n"
  "COMMENT\n"
  "neuromodulation is added as functions:\n"
  "    \n"
  "    modulation = 1 + damod*(maxMod-1)*level\n"
  "\n"
  "where:\n"
  "    \n"
  "    damod  [0]: is a switch for turning modulation on or off {1/0}\n"
  "    maxMod [1]: is the maximum modulation for this specific channel (read from the param file)\n"
  "                e.g. 10% increase would correspond to a factor of 1.1 (100% +10%) {0-inf}\n"
  "    level  [0]: is an additional parameter for scaling modulation. \n"
  "                Can be used simulate non static modulation by gradually changing the value from 0 to 1 {0-1}\n"
  "\n"
  "[] == default values\n"
  "{} == ranges\n"
  "    \n"
  "ENDCOMMENT\n"
  "\n"
  "UNITS {\n"
  "    (mV) = (millivolt)\n"
  "    (mA) = (milliamp)\n"
  "    (S) = (siemens)\n"
  "    (molar) = (1/liter)\n"
  "    (mM) = (millimolar)\n"
  "    FARADAY = (faraday) (coulomb)\n"
  "    R = (k-mole) (joule/degC)\n"
  "}\n"
  "\n"
  "NEURON {\n"
  "    THREADSAFE\n"
  "    SUFFIX cal13\n"
  "    USEION cal READ cali, calo WRITE ical VALENCE 2\n"
  "    RANGE pbar, ical\n"
  "    RANGE damod, maxMod, level, max2, lev2\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    pbar = 0.0 (cm/s)\n"
  "    :q = 1	: room temperature 22-25 C\n"
  "    q = 2	: body temperature 35 C\n"
  "    damod = 0\n"
  "    maxMod = 1\n"
  "    level = 0\n"
  "    max2 = 1\n"
  "    lev2 = 0\n"
  "} \n"
  "\n"
  "ASSIGNED { \n"
  "    v (mV)\n"
  "    ical (mA/cm2)\n"
  "    ecal (mV)\n"
  "    celsius (degC)\n"
  "    cali (mM)\n"
  "    calo (mM)\n"
  "    minf\n"
  "    mtau (ms)\n"
  "    hinf\n"
  "    htau (ms)\n"
  "}\n"
  "\n"
  "STATE { m h }\n"
  "\n"
  "BREAKPOINT {\n"
  "    SOLVE states METHOD cnexp\n"
  "    ical = pbar*m*m*h*ghk(v, cali, calo) *modulation()\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "    rates()\n"
  "    m = minf\n"
  "    h = hinf\n"
  "}\n"
  "\n"
  "DERIVATIVE states { \n"
  "    rates()\n"
  "    m' = (minf-m)/mtau*q\n"
  "    h' = (hinf-h)/htau*q\n"
  "}\n"
  "\n"
  "PROCEDURE rates() {\n"
  "    UNITSOFF\n"
  "    minf = 1/(1+exp((v-(-33))/(-6.7)))\n"
  "    :mtau = 0.06+1/(0.06*exp((v-(-46))/20)+0.41*exp((v-26)/-48))\n"
  "    mtau = 0.06+1/(exp((v-10)/20)+exp((v-(-17))/-48))\n"
  "    hinf = 1/(1+exp((v-(-13.4))/11.9))\n"
  "    htau = 44.3\n"
  "    UNITSON\n"
  "}\n"
  "\n"
  "FUNCTION ghk(v (mV), ci (mM), co (mM)) (.001 coul/cm3) {\n"
  "    LOCAL z, eci, eco\n"
  "    z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))\n"
  "    if(z == 0) {\n"
  "        z = z+1e-6\n"
  "    }\n"
  "    eco = co*(z)/(exp(z)-1)\n"
  "    eci = ci*(-z)/(exp(-z)-1)\n"
  "    ghk = (1e-3)*2*FARADAY*(eci-eco)\n"
  "}\n"
  "\n"
  "FUNCTION modulation() {\n"
  "    : returns modulation factor\n"
  "    \n"
  "    modulation = 1 + damod * ( (maxMod-1)*level + (max2-1)*lev2 ) \n"
  "    if (modulation < 0) {\n"
  "        modulation = 0\n"
  "    }  \n"
  "}\n"
  "\n"
  "COMMENT\n"
  "\n"
  "Activation curve was reconstructed for cultured NAc neurons from\n"
  "P5-P32 Charles River rat pups [1] and shifted to match LVA data [7,\n"
  "Fig.1D]. Activation time constant is from the rodent neuron culture (both\n"
  "rat and mouse cells), room temperature 22-25 C [2, Fig.15A]. Inactivation\n"
  "curve of CaL v1.3 current was taken from HEK cells [3, Fig.2 and p.819]\n"
  "at room temperature.\n"
  "\n"
  "Original NEURON model by Wolf (2005) [4] was modified by Alexander Kozlov\n"
  "<akozlov@kth.se>. Kinetics of m2h type was used [5,6]. Activation\n"
  "time constant was refitted to avoid singularity.\n"
  "\n"
  "[1] Churchill D, Macvicar BA (1998) Biophysical and pharmacological\n"
  "characterization of voltage-dependent Ca2+ channels in neurons isolated\n"
  "from rat nucleus accumbens. J Neurophysiol 79(2):635-47.\n"
  "\n"
  "[2] Kasai H, Neher E (1992) Dihydropyridine-sensitive and\n"
  "omega-conotoxin-sensitive calcium channels in a mammalian\n"
  "neuroblastoma-glioma cell line. J Physiol 448:161-88.\n"
  "\n"
  "[3] Bell DC, Butcher AJ, Berrow NS, Page KM, Brust PF, Nesterova A,\n"
  "Stauderman KA, Seabrook GR, Nurnberg B, Dolphin AC (2001) Biophysical\n"
  "properties, pharmacology, and modulation of human, neuronal L-type\n"
  "(alpha(1D), Ca(V)1.3) voltage-dependent calcium currents. J Neurophysiol\n"
  "85:816-827.\n"
  "\n"
  "[4] Wolf JA, Moyer JT, Lazarewicz MT, Contreras D, Benoit-Marand M,\n"
  "O'Donnell P, Finkel LH (2005) NMDA/AMPA ratio impacts state transitions\n"
  "and entrainment to oscillations in a computational model of the nucleus\n"
  "accumbens medium spiny projection neuron. J Neurosci 25(40):9080-95.\n"
  "\n"
  "[5] Evans RC, Morera-Herreras T, Cui Y, Du K, Sheehan T, Kotaleski JH,\n"
  "Venance L, Blackwell KT (2012) The effects of NMDA subunit composition on\n"
  "calcium influx and spike timing-dependent plasticity in striatal medium\n"
  "spiny neurons. PLoS Comput Biol 8(4):e1002493.\n"
  "\n"
  "[6] Tuckwell HC (2012) Quantitative aspects of L-type Ca2+ currents. Prog\n"
  "Neurobiol 96(1):1-31.\n"
  "\n"
  "[7] Xu W, Lipscombe D (2001) Neuronal cav1.3 L-type channels activate\n"
  "at relatively hyperpolarized membrane potentials and are incompletely\n"
  "inhibited by dihydropyridines. J Neurosci 21(16): 5944-5951.\n"
  "\n"
  "\n"
  "ENDCOMMENT\n"
  ;
#endif
