from kivymd.app import MDApp
from kivy.uix.floatlayout import FloatLayout 
from kivy.garden.matplotlib.backend_kivyagg import FigureCanvasKivyAgg
from skyfield.api import load, load_file, Star  
from skyfield.data import hipparcos, stellarium
from skyfield.magnitudelib import planetary_magnitude
from skyfield.projections import build_stereographic_projection
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

ts = load.timescale()
t = ts.now()

eph = load_file('de421.bsp')
earth = eph['earth']
moon = eph['moon']
sun = eph['sun']
moon_ms = 15**2
sun_ms = 15**2

with load.open('hip_main.dat') as f:
    stars = hipparcos.load_dataframe(f)

with load.open('constellationship.fab') as f:
    constellations = stellarium.parse_constellations(f)

with load.open('star_names.fab') as f:
    star_names = stellarium.parse_star_names(f)


planets = {"mars": "mars barycenter", "mercury": "mercury barycenter", 
           "venus": "venus barycenter", "jupiter": "jupiter barycenter",
           "saturn": "saturn barycenter", "uranus": "uranus barycenter", 
           "neptune": "neptune barycenter"}

all_colors = {"sun": "orange", "moon": "silver", "mercury": "lightgray", 
                 "venus": "yellow", "mars": "darksalmon", "jupiter": "lightsalmon",
                 "saturn": "gold", "uranus": "aqua", "neptune": "deepskyblue"}

edges = [edge for name, edges in constellations for edge in edges]
edges_star1 = [star1 for star1, star2 in edges]
edges_star2 = [star2 for star1, star2 in edges]
star_positions = earth.at(t).observe(Star.from_dataframe(stars))

class MyLayout(FloatLayout):

    plot_exists = False

    def show_planet(self):
        self.ids.error.text = ""
        if self.plot_exists == True:
            self.ids.show_plot.clear_widgets()
            plt.clf()
        try:
            planet_name = self.ids.planet.text
            planet_ = eph[planets[planet_name.lower()]]
            center = earth.at(t).observe(planet_)
            planet_mag = planetary_magnitude(center.apparent())
            projection = build_stereographic_projection(center)
            stars['x'], stars['y'] = projection(star_positions)
            limiting_mag = 6.0
            bright_stars = (stars.magnitude < limiting_mag)
            magnitude = stars['magnitude'][bright_stars]
            ms1 = (0.75 + limiting_mag - magnitude) ** 2
            ms2 = (np.max([6.75 - planet_mag, 0]) + 0.5) ** 2
            xy1 = stars[['x', 'y']].loc[edges_star1].values
            xy2 = stars[['x', 'y']].loc[edges_star2].values
            lines = np.rollaxis(np.array([xy1, xy2]), 1)
            field_of_view = 50.0
            angle = np.pi - field_of_view / 360.0 * np.pi
            limit = np.sin(angle)/(1.0 - np.cos(angle))
            fig, ax = plt.subplots(figsize=(8,8))
            ax.add_collection(LineCollection(lines, color='b'))
            ax.scatter(stars['x'][bright_stars], stars['y'][bright_stars], s=ms1, color='r')
            ax.scatter([0.], [0.], s=ms2, color=all_colors[planet_name.lower()]) 
            ax.annotate(planet_name, (0., 0.))
            for other_planet in planets.keys():
                if other_planet != planet_name.lower():
                    other_planet_ = eph[planets[other_planet]]
                    other_place = earth.at(t).observe(other_planet_)
                    other_x, other_y = projection(other_place)
                    if -limit <= other_x <= limit and -limit <= other_y <= limit:
                        other_mag = planetary_magnitude(other_place.apparent())
                        ms_ = (np.max([6.75 - other_mag, 0]) + 0.5) ** 2
                        ax.scatter([other_x], [other_y], s=ms_, color=all_colors[other_planet])
                        ax.annotate(other_planet[0].upper()+other_planet[1:], (other_x, other_y))
            moon_place = earth.at(t).observe(moon)
            moon_x, moon_y = projection(moon_place)
            if -limit <= moon_x <= limit and -limit <= moon_y <= limit:
                ax.scatter([moon_x], [moon_y], s=moon_ms, color=all_colors['moon'])
                ax.annotate("Moon", (moon_x,moon_y))
            sun_place = earth.at(t).observe(sun)
            sun_x, sun_y = projection(sun_place)
            if -limit <= sun_x <= limit and -limit <= sun_y <= limit:
                ax.scatter([sun_x], [sun_y], s=sun_ms, color=all_colors['sun'])
                ax.annotate("Sun", (sun_x,sun_y))
            nstar_consts = []
            for i in range(87):
                nstar_consts.append(len(constellations[i][1]))
            star_inds = 0
            for i, nstar_const in enumerate(nstar_consts):
                const_position = np.mean(xy1[star_inds:star_inds+nstar_const,:], 0)
                star_inds += nstar_const
                if -limit <= const_position[0] <= limit and -limit <= const_position[1] <= limit:
                    ax.annotate(constellations[i][0], (const_position[0], const_position[1]))
            ax.set_xlim(-limit, limit)
            ax.set_ylim(-limit, limit)
            canvas = FigureCanvasKivyAgg(fig, pos_hint={'x': 0.2, 'y': 0.1}, size_hint=(0.5,0.5))
            self.ids.show_plot.add_widget(canvas)
            canvas.draw()
            self.plot_exists = True
        except Exception:
            self.ids.error.text = "Error!"

    def remove_plot(self):
        self.ids.error.text = ""
        self.ids.show_plot.clear_widgets()
        plt.clf()
        self.plot_exists = False


class MainApp(MDApp):
    def build(self):
        return MyLayout()

app = MainApp()
app.run()