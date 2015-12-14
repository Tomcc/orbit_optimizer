extern crate rand;
extern crate gnuplot;
extern crate cgmath;

use std::cmp;
use std::cmp::Ordering;
use gnuplot::*;
use std::f32;
use std::collections::VecDeque;
use rand::*;
use cgmath::{Vector, Vector2, EuclideanVector};

const GRAVITY:f32 = 9.81;
const COEFFICIENTS:usize = 4;
const SIM_STEP:f32 = 0.1;


fn lerp(src: f32, dst: f32, alpha: f32) -> f32 {
	src + (dst - src) * alpha
}

fn atmospheric_pressure(h: f32) -> f32 {
	let H = 70000.;
	if h >= H {
		0.
	}
	else {
		f32::consts::E.powf(-h/70000.)
	}
}

struct Vessel {
    max_thrust: f32,
    max_burn: f32,
    dry_mass: f32,
    fuel_mass: f32,
    cross_section: f32,
}

impl Vessel {
	fn mass(&self) -> f32 {
		self.dry_mass + self.fuel_mass
	}

	fn find_throttle_for_twr(&self, twr: f32) -> f32 {
		(twr * self.mass() * GRAVITY) / self.max_thrust
	}
}

#[derive(Clone,Debug)]
struct Control {
    coeff: [f32; COEFFICIENTS],
    twr: f32,
    fitness: Option<f32>
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

    fn lerp(a: &Self, b:&Self, s: f32) -> Self {
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

    fn calc_angle(&self, t: f32) -> f32{
    	let mut control_angle = f32::consts::FRAC_PI_2; //start vertical

    	for i in 0..self.coeff.len() {	
    		control_angle += self.coeff[i] * t.powf(i as f32 + 1.);
    	}
    	control_angle
    }

    fn calc_endpoint(&self, vessel: &Vessel) -> Trajectory {

    	let mut trajectory = Trajectory::new(vessel);

    	//advance
    	while trajectory.step(self, vessel, SIM_STEP) {}

    	trajectory
    }

    fn calc_fitness(&mut self, vessel: &Vessel) {
    	let target_y = 71000.;

    	let endpoint = self.calc_endpoint(vessel);

    	let distance_from_target = (target_y - endpoint.pos.y).abs();

    	let f = distance_from_target - endpoint.vel.x;

    	self.fitness = if f.is_normal() && endpoint.pos.y > 10. {
    		Some(f)
    	}
    	else {
    		Some(f32::MAX)
    	}
    }
}

#[derive(Debug)]
struct Trajectory {
	pos: Vector2<f32>,
	vel: Vector2<f32>,
	fuel_mass: f32,
	t: f32,
}

impl Trajectory {
	fn new(vessel: &Vessel) -> Self {
		Trajectory {
			pos: Vector2::new(0.,0.),
			vel: Vector2::new(0.,0.),
			t: 0.,
			fuel_mass: vessel.fuel_mass,
		}
	}

	fn step(&mut self, control: &Control, vessel: &Vessel, dt: f32) -> bool {
		
		self.t += dt;

		let drag = if self.vel.length2() > 0. && false /*HACK*/ {
			let V = self.vel.length();
			let p = atmospheric_pressure(self.pos.y);

			let d = p * V * V * vessel.cross_section * 0.00001;

			-self.vel.normalize() * d
		}
		else {
			Vector2::new(0.,0.)
		};

		let thrust = if self.fuel_mass > 0. {
			let throttle = vessel.find_throttle_for_twr(control.twr);

			self.fuel_mass -= throttle * vessel.max_burn * dt;

			let mass = self.fuel_mass + vessel.dry_mass;

			let control_angle = control.calc_angle(self.t);

			let mut accel = (throttle * vessel.max_thrust) / mass;

			Vector2::new(
				accel * control_angle.cos(),
				accel * control_angle.sin(),
			)
		} 
		else {
			Vector2::new(0.,0.)
		};

		self.vel.x += (thrust.x + drag.x) * dt;
		self.vel.y += (thrust.y + drag.y - GRAVITY) * dt;

		self.pos.x += self.vel.x * dt;
		self.pos.y += self.vel.y * dt;

		self.vel.y > 0.
	}
}

struct Plot {
	figure: Figure,
	trajectories: VecDeque<(Vec<f32>,Vec<f32>)>,
}

impl Plot {
    fn new() -> Self {
    	Plot {
    		figure: Figure::new(),
    		trajectories: VecDeque::new(),
    	}
    } 

    fn add_trajectory(&mut self, vessel: &Vessel, control: &Control) {

		let mut trajectory = Trajectory::new(&vessel);

		let mut x_points:Vec<f32> = vec![];
		let mut y_points:Vec<f32> = vec![];

	   	while trajectory.step(&control, &vessel, SIM_STEP) {
		    x_points.push(trajectory.pos.x);
		    y_points.push(trajectory.pos.y);
		}

		if x_points.len() > 2 {
			if self.trajectories.len() > 5 {
				self.trajectories.pop_front();
			}

			self.trajectories.push_back((
				x_points,
				y_points
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
	};

	let size = 20;

	let mut solutions:Vec<Control> = vec![Control::zero(); size];

	let mut best_fitness = f32::MAX;
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
