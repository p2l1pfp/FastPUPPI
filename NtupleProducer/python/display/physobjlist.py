from PhysicsTools.HeppyCore.utils.deltar import *

def ptsorted(x): return sorted(x, key = lambda o : -o.pt())
def drsorted(x,center): return sorted(x, key = lambda o : deltaR(o.eta(),o.phi(),center[0],center[1]))

class PhysObjList:
    def __init__(self, name, objects, drawers=[], views=None, printer = lambda o : "", sort = ptsorted, modifier = None):
        self.name = name
        self.drawers = [d.clone(name) for d in drawers]
        self.views = views
        self.printer = printer
        self.sort = sort
        if modifier != None:
            self.objects = [] 
            for o in objects:
                oc = o.clone()
                self.objects.append(modifier(oc))
        else:
            self.objects = objects
    def draw(self,view):
        if "all" not in self.views and view not in self.views: 
            return
        for mod in reversed(self.drawers): 
            mod.draw(self.objects)
    def write(self,view,log):
        if "all" not in self.views and view not in self.views: 
            return
        log.write(self.name+"\n")
        ptsum = 0.0
        for o in self.sort(self.objects):
            ptsum += o.pt()
            log.write("      %9.3f  %+5.2f %+5.2f  %s\n" % (o.pt(), o.eta(), o.phi(), self.printer(o)))
        log.write("  TOT %8.2f\n\n" % ptsum)
    def addToLegend(self,view,theLegend):
        if "all" not in self.views and view not in self.views: 
            return
        for mod in reversed(self.drawers): 
            if mod.label and mod.tobjForLegend():
                if mod.label in theLegend.names: continue
                theLegend.AddEntry(mod.tobjForLegend(), mod.label, "P")
                theLegend.names.append(mod.label)
    def writeZoom(self,view,center,size,radius,log):
        if "all" not in self.views and view not in self.views: 
            return
        log.write(self.name+"\n")
        ptsum, once = 0.0, False
        for o in drsorted(self.objects, center):
            dr = deltaR(center[0],center[1],o.eta(),o.phi())
            if dr > radius and ptsum > 0 and not once:
                log.write("      -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n")
                once = True
            if abs(o.eta()-center[0])>size or abs(o.phi()-center[1])>size: 
                continue
            log.write("      %9.3f  %+5.2f %+5.2f [ %+5.2f %+5.2f %.3f ] %s\n" % (o.pt(), o.eta(), o.phi(), o.eta()-center[0], o.phi()-center[1], dr, self.printer(o)))
            if dr < radius: 
                ptsum += o.pt()
        log.write("  TOT %8.2f (for a radius of %.3f)\n\n" % (ptsum,radius))

def read(event,tag,handle,filter=lambda p : True):
    event.getByLabel(tag, handle)
    return [ p for p in handle.product() if filter(p)]

def makeAtCalo(track):
    track.setPtEtaPhiM(track.pt(), track.caloEta(), track.caloPhi(), track.mass())
    return track
