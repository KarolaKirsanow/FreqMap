# -*- coding: utf-8 -*-
"""

@author: kirsanow
"""
"""Plots allele frequency data on a plot from the
Matplotlib basemap toolkit"""


from math import cos, sin, atan2, sqrt, pi
import math
import sys
sys.path.append('C:\WinPython-64bit-2.7.5.3\python-2.7.5.amd64\Lib\site-packages\wx-2.8-msw-unicode')
sys.path.append('C:\WinPython-64bit-2.7.5.3\python-2.7.5.amd64\Lib\site-packages')
import wx
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

import pylab as p
from matplotlib.text import Text
import matplotlib.ticker as ticker
#------------------Set Fonts---------------------------------------------------

path = 'C:\\Windows\\Fonts\\times.ttf' #add path to your system font
prop_1 = font_manager.FontProperties(fname=path)

#==============================================================================
#---------------------Create Draggable Text for Labels-------------------------------------------
#==============================================================================
class DragHandler(object):
    """ A simple class to handle Drag n Drop.

    This is a simple example, which works for Text objects only
    """
    def __init__(self, figure=None) :
        """ Create a new drag handler and connect it to the figure's event system.
        If the figure handler is not given, the current figure is used instead
        """

        if figure is None : figure = p.gcf()
        # simple attibute to store the dragged text object
        self.dragged = None

        # Connect events and callbacks
        figure.canvas.mpl_connect("pick_event", self.on_pick_event)
        figure.canvas.mpl_connect("button_release_event", self.on_release_event)

    def on_pick_event(self, event):
        " Store which text object was picked and were the pick event occurs."

        if isinstance(event.artist, Text):
            self.dragged = event.artist
            self.pick_pos = (event.mouseevent.xdata, event.mouseevent.ydata)
        return True

    def on_release_event(self, event):
        " Update text position and redraw"

        if self.dragged is not None :
            old_pos = self.dragged.get_position()
            new_pos = (old_pos[0] + event.xdata - self.pick_pos[0],
                       old_pos[1] + event.ydata - self.pick_pos[1])
            self.dragged.set_position(new_pos)
            self.dragged = None
            p.draw()
        return True


#==============================================================================
# # -------------- Create Draggable Point class for pie charts----------------------------------
#==============================================================================

class DraggablePoint:
    lock = None #only one can be animated at a time
    def __init__(self, point):
        self.point = point
        self.press = None
        self.background = None

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.point.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease = self.point.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion = self.point.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        if event.inaxes != self.point.axes: return
        if DraggablePoint.lock is not None: return
        contains, attrd = self.point.contains(event)
        if not contains: return
        self.press = (self.point.center), event.xdata, event.ydata
        DraggablePoint.lock = self

        # draw everything but the selected rectangle and store the pixel buffer
        canvas = self.point.figure.canvas
        axes = self.point.axes
        self.point.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.point.axes.bbox)

        # now redraw just the rectangle
        axes.draw_artist(self.point)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        if DraggablePoint.lock is not self:
            return
        if event.inaxes != self.point.axes: return
        self.point.center, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        self.point.center = (self.point.center[0]+dx, self.point.center[1]+dy)

        canvas = self.point.figure.canvas
        axes = self.point.axes
        # restore the background region
        canvas.restore_region(self.background)

        # redraw just the current rectangle
        axes.draw_artist(self.point)

        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        'on release we reset the press data'
        if DraggablePoint.lock is not self:
            return

        self.press = None
        DraggablePoint.lock = None

        # turn off the rect animation property and reset the background
        self.point.set_animated(False)
        self.background = None

        # redraw the full figure
        self.point.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.point.figure.canvas.mpl_disconnect(self.cidpress)
        self.point.figure.canvas.mpl_disconnect(self.cidrelease)
        self.point.figure.canvas.mpl_disconnect(self.cidmotion)

#==============================================================================
# # -------------- Create GUI  and Map Panels----------------------------------
#==============================================================================
class MapPanel(wx.Panel):
    
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)

        ## create sizers
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        grid = wx.GridBagSizer(hgap=5, vgap=5)
        hSizer = wx.BoxSizer(wx.HORIZONTAL)

        ## Create FreqMap button
        self.button =wx.Button(self, label="Create FreqMap")
        self.Bind(wx.EVT_BUTTON, self.OnClick,self.button)
        
        ##Create text inouts for file source and map title
        self.FileLoc = wx.StaticText(self, \
                                     label="Please enter the full path to the .csv")
        grid.Add(self.FileLoc, pos=(1,0))
        self.FileName = wx.TextCtrl\
                        (self,value="D:\Google Drive\datasets\NEW\hg_herc2_rev.csv",\
                         size=(200,-1))
        grid.Add(self.FileName, pos=(2,0))
        
        self.MapNameLbl = wx.StaticText(self, \
                                     label="Please enter the title of the map")
        self.MapName = wx.TextCtrl\
                        (self,value="map name",\
                         size=(200,-1))
        grid.Add(self.MapNameLbl, pos=(3,0))
        grid.Add(self.MapName, pos=(4,0))

        ## Create ComboBox
        self.typeList =  ['Near-Sided Perspective Projection',\
                          'Robinson Projection', \
                 'Orthographic Projection', 'Mercator Projection', \
                 'Gnomonic Projection']
        self.lblhear = wx.StaticText(self, label="Choose map type:")
        grid.Add(self.lblhear, pos=(5,0))
        self.edithear = wx.ComboBox(self, size=(200, -1), \
                                    choices=self.typeList, \
                                    style=wx.CB_DROPDOWN)
        grid.Add(self.edithear, pos=(6,0))
        
        #Create Radio Boxes
        self.radioList = ['Blue Marble', 'Etopo', 'Natural Earth', '19000bce', '13000bce', '6000bce', '4000bce']
        self.rb = wx.RadioBox(self, label="Choose map style:", pos=(20, 210),\
                         choices=self.radioList,  majorDimension=3,
                         style=wx.RA_SPECIFY_ROWS)
        grid.Add(self.rb, pos=(8,0), span=(1,2))


        hSizer.Add(grid, 0, wx.ALL, 5)
        mainSizer.Add(hSizer, 0, wx.ALL, 5)
        mainSizer.Add(self.button, 0, wx.CENTER)
        self.SetSizerAndFit(mainSizer)

    def GetFile(self):
        self.data_source=self.FileName.GetValue().encode('utf-8', 'ignore')
        return self.data_source
        
    def GetTitle(self):
        self.map_title=self.MapName.GetValue().encode('utf-8', 'ignore')
        return self.map_title
        
    def GetType(self):
        self.map_type=self.edithear.GetValue().encode('utf-8', 'ignore')
        return self.map_type
    def GetStyle(self):
        self.map_style=self.rb.GetSelection()
        return self.map_style
        
    def OnClick(self,event):
        panel.GetFile()
        panel.GetTitle()
        panel.GetType()        
        panel.GetStyle()
        panel.draw_map()
        panel.MapType()
        panel.MapStyle()
        panel.draw_pie(self)
        panel.plot_alleles()
        return panel.data_source, panel.map_title, panel.map_type, panel.map_style

#==============================================================================
#  Based on calculation method from http://www.geomidpoint.com/calculation.html
#  This assumes that the lats and longs in the csv are in decimal.
#  To convert to radians from decimal, multiply the decimal * pi/180.
#  Then to convert back from radians to decimal, multiply * 180/pi.
#==============================================================================
    def draw_map(self):
        self.data_file=open(panel.data_source)
        self.lats, self.lons = [], []
        self.labels =[]
        for index, line in enumerate(self.data_file):
            if index > 0:
                self.lats.append(float(line.split(',')[8]))
                self.lons.append(float(line.split(',')[9]))
                self.labels.append(line.split(',')[0])
                self.map_title=panel.map_title
                self.max_lat=max(self.lats)
                self.min_lat=min(self.lats)
                self.max_lon=max(self.lons)
                self.min_lon=min(self.lons)

        self.geolocations=zip(self.lats,self.lons)
        x = 0
        y = 0
        z = 0

        for lat, lon in self.geolocations:
            lat = float(lat*(pi/180))
            lon = float(lon*(pi/180))
            x += cos(lat) * cos(lon)
            y += cos(lat) * sin(lon)
            z += sin(lat)
         
        x = float(x / len(self.geolocations))
        y = float(y / len(self.geolocations))
        z = float(z / len(self.geolocations))
        
        (c_lat, c_lon)=(atan2(z, sqrt(x * x + y * y)),atan2(y, x))
        self.center_lat=c_lat*(180/pi)
        self.center_lon=c_lon*(180/pi)

        ##center point of map determined from center_lat and center_long:
        self.lat_0=self.center_lat; self.lon_0=self.center_lon
        ## altitude of camera (in km):
        self.h = 2000.
        ## resolution = None means don't process the boundary datasets:
 

#==============================================================================
# # --- Select Map Type (from Matplotlib Basemap 1.07 documentation)-----------------
#==============================================================================    
    def MapType(self):
        if panel.map_type=='Near-Sided Perspective Projection':
            self.m1 = Basemap(projection='nsper',satellite_height=self.h*1000.,\
                    lon_0=self.lon_0,lat_0=self.lat_0,resolution='h')    
            
            self.m = Basemap(projection='nsper',satellite_height=self.h*1000.,\
                   lon_0=self.lon_0,lat_0=self.lat_0,resolution='h',\
                llcrnrx=0.,llcrnry=0.,urcrnrx=self.m1.urcrnrx/2.,urcrnry=self.m1.urcrnry/2.)
            self.m.drawcoastlines()
            self.m.drawmapboundary(fill_color='white')

        if panel.map_type=='Gnomonic Projection':    
            self.m1 = Basemap(width=15.e6,height=15.e6,\
                        projection='gnom',lat_0=self.lat_0,lon_0=self.lon_0)
              
            self.m = Basemap(width=15.e6,height=15.e6,\
                        projection='gnom',lat_0=self.lat_0,lon_0=self.lon_0)
            
            self.m.drawcoastlines()
            self.m.drawmapboundary(fill_color='white')
            self.m.drawparallels(np.arange(10,90,20))
            self.m.drawmeridians(np.arange(-180,180,30))     
        
        if panel.map_type=='Mercator Projection':                
            self.m1 = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,\
                    llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='h')            
          
        
            self.m = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,\
                    llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='h')
         
            #self.m.drawcoastlines()
            #self.m.drawparallels(np.arange(-90.,91.,30.))
            #self.m.drawmeridians(np.arange(-180.,181.,60.))
            self.m.drawmapboundary(fill_color='white')             
              
        if panel.map_type=='Robinson Projection':                
            self.m1 = Basemap(projection='robin',lon_0=self.lon_0,resolution='h')    
            self.m = Basemap(projection='robin',lon_0=self.lon_0,resolution='h')
            self.m.drawcoastlines()
            self.m.drawparallels(np.arange(-90.,120.,30.))
            self.m.drawmeridians(np.arange(0.,360.,60.))
            self.m.drawmapboundary(fill_color='white')        
        
        if panel.map_type=='Orthographic Projection':
            self.m1 = Basemap(projection='ortho',lon_0=self.lon_0,lat_0=self.lat_0,resolution='h')    
            self.m = Basemap(projection='ortho',lon_0=self.lon_0,lat_0=self.lat_0,resolution='h')
            self.m.drawcoastlines()
            self.m.drawparallels(np.arange(-90.,120.,30.))
            self.m.drawmeridians(np.arange(0.,420.,60.))
            self.m.drawmapboundary(fill_color='white')
            
  ##plot mapshading onto m1, not m!(also centers the map in the window!)           
    def MapStyle(self):
        if panel.map_style==0:
            self.m1.bluemarble()            
        if panel.map_style==1:
            self.m1.etopo()
        if panel.map_style==2:
            self.m1.warpimage(image="D:\Google Drive\BASEmaps\Natural_Earth_hypso\HYP_HR_SR_OB_DR copy.png")#Nat_Earth_hypso 70mb
           #self.m1.shadedrelief()
        if panel.map_style==3:
            self.m1.warpimage(image="D:\\Google Drive\\BASEmaps\\Hypso_19000bce.png")
        if panel.map_style==4:
            self.m1.warpimage(image="D:\\Google Drive\\BASEmaps\\Hypso_13000bce.png")
        if panel.map_style==5:
            self.m1.warpimage(image="D:\\Google Drive\\BASEmaps\\Hypso_6000bce.png")
        if panel.map_style==6:
            self.m1.warpimage(image="D:\\Google Drive\\BASEmaps\\Hypso_4000bce.png")
        else:
            None    
        ## plot site markers:
        self.x,self.y = self.m(self.lons, self.lats)# NB lons precede lats
        ## plot geographical center:
        a,b=self.m(self.center_lon, self.center_lat)
        ## plot map title:
        plt.title(panel.map_title,fontsize=12,family='Times New Roman', weight='demibold')
        ## offset labels from points for legibility:
        for label, xpt, ypt in zip(self.labels, self.x, self.y):
            if panel.map_style==1:
                p.text(xpt+80000, ypt+80000, label, fontsize=14, family='Times New Roman',\
                weight='semibold', color='#463E3F', picker=True,\
                bbox = dict(boxstyle = 'round,pad=0.2', fc='Cornsilk', \
                ec="none", alpha = 0.8)) 
            elif panel.map_style==2:
                p.text(xpt+80000, ypt+80000, label, fontsize=14,family='Times New Roman',\
                weight='semibold', color='#463E3F',picker=True,\
                bbox = dict(boxstyle = 'round,pad=0.2', fc='Cornsilk', \
                ec="none", alpha = 0.8))
            else:
                p.text(xpt+80000, ypt+80000, label,fontsize=14, family='Times New Roman', \
                weight='semibold', color='black',picker=True,\
                bbox = dict(boxstyle = 'round,pad=0.2', fc='Cornsilk', \
                ec="none", alpha = 0.8))


#==============================================================================
# -----------------Draw Allele Frequency Piecharts ----------------------------
# (Modified from Thomas Lecocq @Geophysique.be + Manuel Metz @matplotlib.org)
#==============================================================================
# # The draw_pie method allows a custom number of parts to be set.
# # Each part of the pie is drawn as a scatter point w/ a custom wedge marker.
# # The ratios list is provided to the draw_pie method.
# # X and Y are the coordinates of the center of the pie. 
# # Size is the size of the pie chart. 
# # The first part is drawn from 0.0 to 2*pi*ratios[0], 
# # the second between 2*pi*ratios[0] to 2*pi*ratios[0]+2*pi*ratios[1],...
# # until 2*pi is reached. 
# # The colors of the parts are defined in the colors property 
# # The linspace function generates linearly spaced vectors - this gives direct 
# # control over the number of points. 
# # y = linspace(a,b,n) generates a row vector y of n points linearly spaced
# # between and including a and b. For n < 2, linspace returns b.
#==============================================================================                
    def draw_pie(self,ax, ratios=[0.4,0.3,0.3], X=0, Y=0, size = 0):
        self.colors = ['#F75D59','honeydew','#4E387E','honeydew']
        xy = []
        start = 0.  
        for freq in ratios:
            x = [0] + np.cos(np.linspace(2*math.pi*start,2*math.pi*(start+freq), 30)).tolist()
            y = [0] + np.sin(np.linspace(2*math.pi*start,2*math.pi*(start+freq), 30)).tolist()
            xy.append(zip(x,y))
            start =start+ freq
            #print start
                         
        for i, xyi in enumerate(xy):
            self.ax.scatter([X],[Y], marker=(xyi,0), s=size,  facecolor=self.colors[i],edgecolor='#463E3F', linewidth=2)
    
                 
    fig = plt.figure(figsize=(11.7,8.3))
    ##plot subplots - one row, one column, first plot
    ax = plt.subplot(111)
    
    
    def plot_alleles(self):            
        self.data = np.loadtxt(panel.data_source,skiprows=1,
         dtype={'names': ('group', 'allele1','allele1_freq',\
         'allele2','allele2_freq','allele3','allele3_freq',\
         'gene','lat','lon','no2N'),
         'formats': ('S10', 'S2','float','S2','float','S4','float',\
         'S7','float','float','float')}, delimiter=',')
        
               
        if self.data.ndim==0:
            
            
            X, Y = self.m(self.data['lon'],self.data['lat'])
            
            
            a1=self.data['allele1']
            a2=self.data['allele2']
            a3=self.data['allele3']
            self.r1=self.data['allele1_freq']
            self.r2=self.data['allele2_freq']
            self.r3=self.data['allele3_freq']
            legend_title=self.data['gene']
           
            individuals=self.data['no2N']
            if individuals==1:
                self.draw_pie(ax=plt.subplot(111), ratios=[self.r1,self.r2,self.r3], X=X,Y=Y,size=100)
            else:
                self.draw_pie(ax=plt.subplot(111), ratios=[self.r1,self.r2,self.r3], X=X,Y=Y,size=300)
       
        else:
            
            for di in self.data:
                group, allele1,allele1_freq ,allele2,\
                allele2_freq,allele3,\
                allele3_freq, gene,lat,lon,no2N =di
                   
                X, Y = self.m(lon,lat)
                a1=allele1
                a2=allele2
                a3=allele3
                self.r1=allele1_freq
                self.r2=allele2_freq
                self.r3=allele3_freq
                legend_title=gene
                
                individuals=no2N
                if individuals==1:
                   self.draw_pie(ax=plt.subplot(111), ratios=[self.r1,self.r2,self.r3], X=X,Y=Y,size=no2N*60)
                else:
                   self.draw_pie(ax=plt.subplot(111), ratios=[self.r1,self.r2,self.r3], X=X,Y=Y,size=no2N*60)
        
         ##Create Draggable Pie Charts:
#        for i in group:
#            pies =panel.plot_alleles()
#            drs=[]    
#                    
#            for pie in pies:    
#                dr = DraggablePoint(pie)
#                dr.connect()
#                drs.append(dr)
#                print drs
        
        
        ## handle differences in allele number          
        if a3=="none":
            self.alleles=[a1,a2]
        else:
            self.alleles=[a1,a2,a3]
        ## create legend    
        self.proxies=[]
        for i,ai in enumerate (self.alleles):#need a list to iterate over
            self.proxies.append(plt.Rectangle((0, 0), 1, 1, facecolor=self.colors[i]))
            
        self.ax.legend(self.proxies,self.alleles,title=legend_title)
        print self.map_style, 'ms', self.map_type, 'mt'
        plt.show()            

plt.rc('font',**{'family':'serif','serif':['Times New Roman'],'size':12})
## Create the event handler for Draggable Text:
dragh = DragHandler()

p.show()

if __name__ == '__main__':
    app = wx.App(False)
    frame = wx.Frame(None, title='FreqMap', size =(550,400))
    panel=MapPanel(frame)

    frame.Show()

app.MainLoop()
#robinson projectioon will not plot basemap background - is the input file 
#format somehow the problem?