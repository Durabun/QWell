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


class QwellMenuLayout(GridLayout):
	pass

class QwellMenu(App):
	
	def build(self):
		return QwellMenuLayout()


if __name__ == '__main__':
    QwellMenu().run()
