extern crate rand;
extern crate gnuplot;
extern crate cgmath;

use gnuplot::*;
use std::f64;
use std::collections::VecDeque;
use rand::*;
use cgmath::{Vector, Vector2, EuclideanVector};

const COEFFICIENTS:usize = 2;
const SIM_STEP:f64 = 0.1;
const G:f64 = 6.67384E-11;
const PLANET_MASS:f64 = 5.2915793E22;
const PLANET_RADIUS:f64 = 600000.;
const PLANET_CENTER:Vector2<f64> = Vector2{x: 0., y: PLANET_RADIUS};

fn gravity_at(pos: Vector2<f64>) -> Vector2<f64> {
	Vector2{
		x: 0.,
		y: -((G * PLANET_MASS) / (PLANET_RADIUS + pos.y).powf(2.))
	}
}

fn altitude(pos: Vector2<f64>) -> f64 {
	pos.y
}

fn lerp(src: f64, dst: f64, alpha: f64) -> f64 {
	src + (dst - src) * alpha
}

fn atmospheric_pressure(pos: Vector2<f64>) -> f64 {
	let H = 70000.;
	if pos.y >= H {
		0.
	}
	else {
		f64::consts::E.powf(-pos.y/70000.)
	}
}

struct Vessel {
    max_thrust: f64,
    max_burn: f64,
    dry_mass: f64,
    fuel_mass: f64,
    cross_section: f64,
    starting_height: f64,
}

impl Vessel {
	fn mass(&self) -> f64 {
		self.dry_mass + self.fuel_mass
	}

	fn find_throttle_for_twr(&self, twr: f64, pos: Vector2<f64>) -> f64 {
		assert!(twr > 0.);
		let t = (twr * self.mass() * gravity_at(pos).length()) / self.max_thrust;
		assert!(t <= 1.);
		t
	}
}

#[derive(Clone,Debug)]
struct Control {
    coeff: [f64; COEFFICIENTS],
    twr: f64,
    fitness: Option<f64>
}

impl Control {
    fn zero() -> Self {
    	Control {
    		coeff:[0.;COEFFICIENTS],
    		twr: 1.5, //maybe it should vary?
    		fitness: None
    	}
    }

   	fn duplicate(&self) -> Self {
		Control {
			coeff: self.coeff,
			twr: self.twr,
			fitness: None
		}
	}

    fn lerp(a: &Self, b:&Self, s: f64) -> Self {
    	let mut res = Control::zero();

    	for i in 0..res.coeff.len() {
    		res.coeff[i] = lerp(a.coeff[i], b.coeff[i], s);
    	}

    	res
    }

    fn mutated(&self) -> Self {
    	let mut res = Control::duplicate(self);
    	let ptr: *mut Control = &mut res;
		let bytes: *mut u8 = ptr as *mut u8;
		unsafe {
			let bitlen = std::mem::size_of_val(&self.coeff) * 8;
			let target_bit:usize = thread_rng().gen_range(0, bitlen);
			let byte = target_bit / 8;
			let bit = target_bit % 8;

			(*(bytes.offset(byte as isize))) ^= 1 << bit;
		}

		res
    }

    fn mate(a: &Self, b:&Self) -> (Self, Self) {
    	( 
    		Self::lerp(a, b, 0.25).mutated(), 
    		Self::lerp(b, a, 0.25).mutated() 
    	)
    }

    fn calc_angle(&self, t: f64) -> f64{
    	let mut control_angle = f64::consts::FRAC_PI_2; //start vertical

    	for i in 0..self.coeff.len() {	
    		control_angle += self.coeff[i] * t.powf(i as f64 + 1.);
    	}
    	control_angle
    }

    fn calc_endpoint(&self, vessel: &Vessel) -> Trajectory {

    	let mut trajectory = Trajectory::new(vessel, self);

    	//advance
    	while trajectory.step(self, vessel, SIM_STEP) {}

    	trajectory
    }

    fn calc_fitness(&mut self, vessel: &Vessel) {
    	let target_alt = 71000.;

    	let endpoint = self.calc_endpoint(vessel);

    	let distance_from_target = (target_alt - endpoint.altitude()).abs();

    	let f = distance_from_target - endpoint.tangent_speed();

    	self.fitness = if f.is_normal() && endpoint.altitude() > 10. {
    		Some(f)
    	}
    	else {
    		Some(f64::MAX)
    	}
    }
}

#[derive(Debug)]
struct Trajectory {
	pos: Vector2<f64>,
	vel: Vector2<f64>,
	fuel_mass: f64,
	t: f64,
	throttle: f64,
	control_angle:f64,
}

impl Trajectory {
	fn new(vessel: &Vessel, control: &Control) -> Self {
		let start_pos = Vector2::new(0.,vessel.starting_height);
		Trajectory {
			pos: start_pos,
			vel: Vector2::new(0.,0.),
			t: 0.,
			fuel_mass: vessel.fuel_mass,
			throttle: vessel.find_throttle_for_twr(control.twr, start_pos),
			control_angle: f64::consts::FRAC_PI_2,
		}
	}

	fn altitude(&self) -> f64 {
		altitude(self.pos)
	}

	fn vertical_speed(&self) -> f64 {
		self.vel.y
	}

	fn tangent_speed(&self) -> f64 {
		self.vel.x
	}

	fn step(&mut self, control: &Control, vessel: &Vessel, dt: f64) -> bool {
		
		self.t += dt;

		let drag = if self.vel.length2() > 0. && false /*HACK*/ {
			let V = self.vel.length();
			let p = atmospheric_pressure(self.pos);

			let d = p * V * V * vessel.cross_section * 0.00001;

			-self.vel.normalize() * d
		}
		else {
			Vector2::new(0.,0.)
		};

		let thrust = if self.fuel_mass > 0. {
			self.fuel_mass -= self.throttle * vessel.max_burn * dt;

			let mass = self.fuel_mass + vessel.dry_mass;

			self.control_angle = control.calc_angle(self.t);

			let mut accel = (self.throttle * vessel.max_thrust) / mass;

			Vector2::new(
				accel * self.control_angle.cos(),
				accel * self.control_angle.sin(),
			)
		} 
		else {
			Vector2::new(0.,0.)
		};

		let g = gravity_at(self.pos);

		self.vel.x += (thrust.x + drag.x + g.x) * dt;
		self.vel.y += (thrust.y + drag.y + g.y) * dt;

		self.pos.x += self.vel.x * dt;
		self.pos.y += self.vel.y * dt;

		self.vertical_speed() > 0.
	}
}

struct Plot {
	figure: Figure,
	trajectories: VecDeque<(Vec<f32>,Vec<f32>,Vec<f32>,Vec<f32>)>,
}

impl Plot {
    fn new() -> Self {
    	Plot {
    		figure: Figure::new(),
    		trajectories: VecDeque::new(),
    	}
    } 

    fn add_trajectory(&mut self, vessel: &Vessel, control: &Control) {

    	let RADDEG = 180. / f64::consts::PI;

		let mut trajectory = Trajectory::new(&vessel, &control);

		let mut x_points:Vec<f32> = vec![];
		let mut y_points:Vec<f32> = vec![];
		let mut a_points:Vec<f32> = vec![];
		let mut f_points:Vec<f32> = vec![];

	   	while trajectory.step(&control, &vessel, SIM_STEP) {
		    x_points.push(trajectory.pos.x as f32);
		    y_points.push(trajectory.pos.y as f32);
		    a_points.push((trajectory.control_angle * RADDEG * 1000.) as f32);
			f_points.push(trajectory.fuel_mass as f32 * 10.);
		}

		if x_points.len() > 2 {
			if self.trajectories.len() > 5 {
				self.trajectories.pop_front();
			}

			self.trajectories.push_back((
				x_points,
				y_points,
				a_points,
				f_points,
			));
		}
    }


	fn refresh(&mut self) {
		self.figure.clear_axes();
		{
			let mut plot = self.figure.axes2d();

			let r = 11./16.;
			plot.set_aspect_ratio(Fix(r as f64));
	    	plot.set_x_range(Fix(-100.0), Fix(140000.0));
	    	plot.set_y_range(Fix(-100.0), Fix(140000.0 * r));

			for t in &self.trajectories {
				plot.lines(
					&t.0,
					&t.1,
					&[
						Caption("Best ascent"),
						Color("black")
					]
				);

				plot.lines(
					&t.0,
					&t.2,
					&[
						Caption("Angles"),
						Color("red")
					]
				);

				plot.lines(
					&t.0,
					&t.3,
					&[
						Caption("Fuel"),
						Color("blue")
					]
				);
			}
		}

		self.figure.show();
	}
}

fn sort_solutions(solutions: &mut [Control], vessel: &Vessel) {
	for i in 0..solutions.len() {
		solutions[i].calc_fitness(&vessel);
	}

	solutions.sort_by(|a, b| {
		a.fitness.unwrap().partial_cmp(&b.fitness.unwrap()).unwrap()
	});
}

fn main() {
	let mut plot = Plot::new();

	let vessel = Vessel {
		max_thrust: 170182.54,
		max_burn: 63.98,
		dry_mass: 2776.32,
		fuel_mass: 4035.,
		cross_section: 1.5,
		starting_height: 0.,
	};

	let size = 20;

	let mut solutions:Vec<Control> = vec![Control::zero(); size];

	let mut best_fitness = f64::MAX;
	loop {
		sort_solutions(&mut solutions, &vessel);

		if let Some(fitness) = solutions[0].fitness {
			if fitness < best_fitness {
				best_fitness = fitness;
				println!("New best fitness: {}", best_fitness);
				println!("{:?}", solutions[0]);

				plot.add_trajectory(&vessel, &solutions[0]);
				plot.refresh();
			}
		}

		//reproduce the first half

		let mut new_gen = Vec::with_capacity(size);
			for i in 0..size/2 {
			let j:usize = thread_rng().gen_range(0, size/2);

			let spawn = Control::mate(
				&solutions[i],
				&solutions[j]
			);

			new_gen.push(spawn.0);
			new_gen.push(spawn.1);
		}

		solutions = new_gen;
	}
}
