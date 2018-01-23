import sys
import json

d = json.load(open('data/source.json'))
nobs   = d['bindata']['data'][0]
b      = d['bindata']['bkg'][0]
deltab = d['bindata']['bkgerr'][0]
s      = d['bindata']['sig'][0]

# derived data
tau = b/deltab/deltab
mobs = round(tau*b)

print 'tau: ', tau, 'm: ', mobs

import ROOT
import sys
w = ROOT.RooWorkspace("w",True);

#-----------------

w.factory("prod:nsig(mu[1,0,10],s[1])");
w.factory("sum:nexp_sr(nsig,b[1,40,300])");
w.factory("Poisson:on_model(nobs_sr[0,1000],nexp_sr)");

#-----------------

w.var('s').setVal(s)
w.var('b').setVal(b)

w.var('s').setConstant(True)
w.var('nobs_sr').setVal(nobs)


w.factory("prod:nexp_cr(tau[1],b)");
w.factory("Poisson:off_model(nobs_cr[0,1000],nexp_cr)");
w.var('nobs_cr').setVal(mobs)
w.var('nobs_cr').setConstant(True)
w.var('tau').setVal(tau)
w.var('tau').setConstant(True)

w.factory("PROD:onoff(on_model,off_model)");


data = ROOT.RooDataSet('data','data', ROOT.RooArgSet(w.var('nobs_sr'), w.var('nobs_cr')))
data.add(ROOT.RooArgSet(w.var('nobs_sr'), w.var('nobs_cr')))

getattr(w,'import')(data)






modelConfig = ROOT.RooStats.ModelConfig(w);
modelConfig.SetPdf(w.pdf('onoff'));
modelConfig.SetParametersOfInterest(ROOT.RooArgSet(w.var('mu')))
modelConfig.SetNuisanceParameters(ROOT.RooArgSet(w.var('b')))
modelConfig.SetObservables(ROOT.RooArgSet(w.var('nobs_sr'), w.var('nobs_cr')))
modelConfig.SetGlobalObservables( ROOT.RooArgSet())
modelConfig.SetName("ModelConfig")
getattr(w,'import')(modelConfig)

w.Print()

for v in ('b','mu','nobs_cr','nobs_sr','s','tau'):
    w.var(v).Print()

nll = w.pdf('onoff').createNLL(data)

print 'nll', nll.getVal()
print 'lambda', 2*nll.getVal()

# print 'ok'
data.get(0).Print('nobs_sr nobs_cr')


sbModel = w.obj('ModelConfig')
poi = sbModel.GetParametersOfInterest().first()

sbModel.SetSnapshot(ROOT.RooArgSet(poi))

bModel = sbModel.Clone()
bModel.SetName("bonly")
poi.setVal(0)
bModel.SetSnapshot(ROOT.RooArgSet(poi))


ac = ROOT.RooStats.AsymptoticCalculator(data, bModel, sbModel)
ac.SetOneSided(True)
# ac.SetQTilde(False)
ac.SetPrintLevel(10)
ROOT.RooStats.AsymptoticCalculator.SetPrintLevel(10)


print 'asimov'
poi.setVal(0.0)
poi.setConstant(True)

poialt = ac.GetAlternateModel().GetNuisanceParameters()
tmp  = poialt.snapshot()
roofit_asimov = ac.MakeAsimovData(data, ac.GetNullModel(), poialt,tmp)
print 'roofit asimov', roofit_asimov.Print('v')


manual_asimov = ROOT.RooDataSet('data','data', ROOT.RooArgSet(w.var('nobs_sr'), w.var('nobs_cr')))
basimov_manual = 8.27528e+01
nasimov_manual = basimov_manual
masimov_manual = basimov_manual*tau

w.var('nobs_sr').setVal(nasimov_manual)
w.var('nobs_cr').setVal(masimov_manual)
manual_asimov.add(ROOT.RooArgSet(w.var('nobs_sr'), w.var('nobs_cr')))

print 'manual asimov: ', manual_asimov.Print('v')


print '-----------'

poi.setVal(2.0)

   #
   # 1  b            8.00153e+01   9.93330e+00   2.13054e-04   3.46241e-02
   # 2  mu           5.02371e-01   1.17946e+00   1.13016e-03   6.89767e-03

  # NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
  #  1  b            8.27526e+01   6.76409e+00   2.08970e-04  -3.14627e-02


nll = w.pdf('onoff').createNLL(data)
asimov_data = manual_asimov


asimov_nll = w.pdf('onoff').createNLL(asimov_data)


bestfitb  = 8.00153e+01
bestfitmu = 5.02371e-01

poi.setVal(bestfitmu)
w.var('b').setVal(bestfitb)

bestfitnll_observed = nll.getVal()
bestfitnll_asimov = asimov_nll.getVal()

condfitb  = 7.25758e+01
condfitmu = 2.0

poi.setVal(condfitmu)
w.var('b').setVal(condfitb)

condfitnll_observed = nll.getVal()
condfitnll_asimov = asimov_nll.getVal()


print 'mytau:', w.var('tau').getVal()
print 'mynll ASIMOV',   bestfitnll_asimov, condfitnll_asimov
print 'mynll OBSERVED', bestfitnll_observed, condfitnll_observed

print '--- ASIMOV FIT --'
poi.setVal(0)
poi.setConstant(True)
w.pdf('onoff').fitTo(asimov_data)

print 'ASHA_NLL', w.pdf('onoff').createNLL(asimov_data).getVal()
print 'ASHA_WORKS', asimov_nll.getVal()
print 'ASHA', w.var('b').getVal()

print '------'
poi.setVal(0.4)
poi.setConstant(False)
sbModel.SetSnapshot(ROOT.RooArgSet(poi))
r = ac.GetHypoTest()
r.SetBackgroundAsAlt(True)
r.Print()

import math
from scipy.stats import norm
for nsigma in [-2,-1,0,1,2]:
    print nsigma,'null', r.NullPValue()
    print nsigma,'alt',  r.AlternatePValue()
    qmu = ROOT.Math.normal_quantile_c( r.NullPValue(),1.)**2
    qmuA = (ROOT.Math.normal_quantile( r.AlternatePValue(),1.) + math.sqrt(qmu))**2
    print nsigma,'qmu', qmu
    print nsigma,'qmuA', qmuA

    qmu = (math.sqrt(qmuA)-nsigma)**2
    CLsb = 1-norm.cdf(math.sqrt(qmu) )
    CLb  = norm.cdf(math.sqrt(qmuA)-math.sqrt(qmu))
    print nsigma,'exp', CLsb/CLb

    # print ac.GetExpectedPValues( r.AlternatePValue(), r.NullPValue(), nsigma, True, True);
    print nsigma,ac.GetExpectedPValues( r.NullPValue(), r.AlternatePValue(), nsigma, True, True);

print 'hmm', r.CLs()
# sys.exit(0)

print 'muhat'
ac.GetMuHat().Print()

poi.setVal(0.4)

# sys.exit(0)

h = nll.createHistogram('mu,b',50,50)
# fr = ROOT.RooPlot(w.var('mu'), w.var('b'))
# nll.plotOn(fr)
c = ROOT.TCanvas('c','c',800,800)
h.Draw('cont4z')
# fr.Draw()
c.SaveAs('c.pdf')


prof = nll.createProfile(ROOT.RooArgSet(w.var('mu')))
fr = w.var('mu').frame(0,2)
nll.plotOn(fr,ROOT.RooFit.ShiftToZero())
prof.plotOn(fr, ROOT.RooFit.LineColor(ROOT.kRed))
fr.Draw()
c.SaveAs('p.pdf')

w.var('b').setVal(4.99777e+01)
w.var('mu').setVal(1.00226e+00)
print nll.getVal()


nll = w.pdf('onoff').createNLL(data)
w.var('b').setVal(5.49355e+01)
w.var('mu').setVal(0.4)
print nll.getVal()

print 'done'
# ###
import sys
# sys.exit(0)





calc = ROOT.RooStats.HypoTestInverter(ac)
calc.RunFixedScan(51,0,5)
calc.SetConfidenceLevel(0.95)
calc.UseCLs(True)

result = calc.GetInterval()

plot = ROOT.RooStats.HypoTestInverterPlot("plot","plot",result)
c = ROOT.TCanvas()
c.SetLogy(False)
plot.Draw("OBS EXP CLb 2CL")
c.Draw()
c.SaveAs('scan.pdf')

print 'observed: ', result.UpperLimit()
for i in [-2,-1,0,1,2]:
    print 'expected {}: '.format(i), result.GetExpectedUpperLimit(i)
