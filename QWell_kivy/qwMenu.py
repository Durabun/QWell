import kivy
#kivy.require('1.0.6') # replace with your current kivy version !

from kivy.app import App
from kivy.uix.label import Label
from kivy.uix.gridlayout import GridLayout
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.button import Button
from kivy.uix.widget import Widget
from kivy.uix.textinput import TextInput
from kivy.uix.anchorlayout import AnchorLayout
from kivy.uix.popup import Popup
import matplotlib.pyplot as plt
import numpy as np


class QwellMenu(App):
	
	def build(self):
		
		layout = GridLayout(cols = 2)
		layout.add_widget(Label(text = 'Mass'))
		self.Mass = TextInput(multiline = False)
		layout.add_widget(self.Mass)
		layout.add_widget(Label(text = 'Depth'))
		self.Depth = TextInput(multiline = False)
		layout.add_widget(self.Depth)
		layout.add_widget(Label(text = 'Width'))
		self.hWidth = TextInput(multiline = False)
		layout.add_widget(self.hWidth)
		self.popup = Popup(title = 'Test popup')
		GoBtn = Button(text = 'GO!')
#		GoBtn.bind(on_release = self.callback)
		GoBtn.bind(on_release = self.plot)
		layout.add_widget(GoBtn)
		anchor = AnchorLayout(anchor_x = 'center', anchor_y = 'bottom')
		btn = Button(text = ('test'))
		anchor.add_widget(btn)	
		

		return layout
		return anchor
	
	def callback(self,btn):
		print("Mass is "+self.Mass.text)
		print("Depth is "+self.Depth.text)
		print("Half Width is "+self.hWidth.text)
	
	def plot(self,btn):
		m = float(self.Mass.text)
		d = float(self.Depth.text)
		a = float(self.hWidth.text)
		t = np.arange(0.0,2.0,0.01)
		s = a+np.sin(m*d*np.pi*t)
		plt.plot(t,s)
		return plt.show()


if __name__ == '__main__':
    QwellMenu().run()
