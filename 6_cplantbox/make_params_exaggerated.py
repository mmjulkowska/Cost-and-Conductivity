"""Exaggerated topologies (requested by M. Julkowska, 20 Nov 2026): fewer and longer broomstick laterals, a stronger length gradient
in the Christmas tree, uniform lateral length in the droopy telephone pole. Hand-set from the v1 fit (values found with try_exag.py);
lateral angles unchanged. Writes params_exaggerated.json."""
import cpb_paths as CP   # paths; changes into outputs/cplantbox/
import json
F=json.load(open('params_fitted.json')); E={}
ov={'B':dict(r=1.55,la=4.0,ln=0.25+0.1*(15-4.0-0.1)/2,lnk=0.1,rL=3.0,lmaxL=10.),
    'C':dict(r=2.6,la=1.0,ln=0.5,lnk=0.,rL=1.4,lmaxL=10.),
    'dT':dict(r=2.2,la=0.8,ln=0.38,lnk=0.,rL=3.,lmaxL=1.2,gfL=2)}
for t in F:
    E[t]={g:{**F[t][g],**ov[t]} for g in F[t]}
json.dump(E,open('params_exaggerated.json','w'),indent=1)
