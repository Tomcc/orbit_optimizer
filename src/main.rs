extern crate rand;
extern crate gnuplot;

use std::cmp::Ordering;
use gnuplot::*;
use std::f32;
use std::collections::VecDeque;

const GRAVITY:f32 = 9.81;


fn lerp(src: f32, dst: f32, alpha: f32) -> f32 {
	src + (dst - src) * alpha
}

fn random() -> f32 {
	(rand::random::<f32>() * 2.) - 1.
}

struct Vessel {
    max_thrust: f32,
    max_burn: f32,
    dry_mass: f32,
    fuel_mass: f32
}

impl Vessel {
	fn mass(&self) -> f32 {
		self.dry_mass + self.fuel_mass
	}

	fn find_throttle_for_TWR(&self, twr: f32) -> f32 {
		(twr * self.mass() * GRAVITY) / self.max_thrust
	}
}

#[derive(Clone, Debug)]
struct Control {
    coeff: [f32; 2],
    twr: f32,
    fitness: Option<f32>
}

impl Control {
    fn zero() -> Self {
    	Control {
    		coeff:[0.;2],
    		twr: 1.5, //maybe it should vary?
    		fitness: None
    	}
    }

    fn random() -> Self {
    	let mut s = Control::zero();
    	for c in &mut s.coeff {
    		*c = random();
    	}
    	s
    }

    fn lerp(a: &Self, b:&Self, s: f32) -> Self {
    	let mut res = Control::zero();

    	for i in 0..res.coeff.len() {
    		res.coeff[i] = lerp(a.coeff[i], b.coeff[i], s);
    	}

    	res
    }

    fn mutated(&self) -> Self {
    	let mut res = Control::zero();

    	for i in 0..res.coeff.len() {
    		res.coeff[i] = self.coeff[i] + random() * 0.001
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
    	let mut control_angle = f32::consts::PI / 2.; //start vertical

    	for i in 0..self.coeff.len() {
    		control_angle += self.coeff[i] * t.powf(i as f32);
    	}
    	control_angle
    }

    fn calc_endpoint(&self, vessel: &Vessel) -> Trajectory {

    	let mut trajectory = Trajectory::new(vessel);

    	//advance
    	while trajectory.step(self, vessel, 1.0) {}

    	trajectory
    }

    fn calc_fitness(&mut self, vessel: &Vessel) {
    	let target_y = 71000.;

    	let endpoint = self.calc_endpoint(vessel);

    	let distance_from_target = (target_y - endpoint.y).abs();

    	self.fitness = Some(distance_from_target - endpoint.v_x);
    }

    fn fitness(&self) -> f32 {
    	self.fitness.unwrap()
    }
}


struct Trajectory {
	y: f32,
	x: f32,
	v_x: f32,
	v_y: f32,
	fuel_mass: f32,
	t: f32,
}

impl Trajectory {
	fn new(vessel: &Vessel) -> Self {
		Trajectory {
			x: 0.,
			y: 0.,
			v_x: 0.,
			v_y: 0.,
			t: 0.,
			fuel_mass: vessel.fuel_mass,
		}
	}

	fn step(&mut self, control: &Control, vessel: &Vessel, dt: f32) -> bool {
		let throttle = if self.fuel_mass > 0. {
			vessel.find_throttle_for_TWR(control.twr)
		} 
		else {
			0.
		};


		self.t += dt;

		self.fuel_mass -= throttle * vessel.max_burn * dt;

		let mass = self.fuel_mass + vessel.dry_mass;
		let control_angle = control.calc_angle(self.t);

		let accel = (throttle * vessel.max_thrust * dt) / mass;

		self.v_x += accel * control_angle.cos() * dt;
		self.v_y += (accel * control_angle.sin() - GRAVITY) * dt;

		self.x += self.v_x * dt;
		self.y += self.v_y * dt;

		self.v_y > 0.
	}
}

struct Plot {
	figure: Figure,
	trajectories: VecDeque<(Vec<f32>,Vec<f32>)>,
}

impl Plot {
    fn new() -> Self {
    	let mut fg = Figure::new();
    	Plot {
    		figure: Figure::new(),
    		trajectories: VecDeque::new(),
    	}
    } 

    fn add_trajectory(&mut self, vessel: &Vessel, control: &Control) {

		let mut trajectory = Trajectory::new(&vessel);

		let mut x_points:Vec<f32> = vec![];
		let mut y_points:Vec<f32> = vec![];

	   	while trajectory.step(&control, &vessel, 1.0) {
		    x_points.push(trajectory.x);
		    y_points.push(trajectory.y);
		}

		if self.trajectories.len() > 5 {
			self.trajectories.pop_front();
		}

		self.trajectories.push_back((
			x_points,
			y_points
		));
    }


	fn refresh(&mut self) {
		self.figure.clear_axes();
		{
			let mut plot = self.figure.axes2d();
	    	plot.set_x_range(Fix(0.0), Fix(160000.0));
	    	plot.set_y_range(Fix(0.0), Fix(160000.0));

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
		a.fitness().partial_cmp(&b.fitness()).unwrap()
	});

	assert!(solutions[0].fitness() <= solutions[solutions.len()-1].fitness());
}

fn main() {
	let mut plot = Plot::new();

	let vessel = Vessel {
		max_thrust: 170182.54,
		max_burn: 63.98,
		dry_mass: 2776.32,
		fuel_mass: 4035.,
	};

	let size = 500;

	let mut solutions:Vec<Control> = vec![Control::zero(); size];
	// for s in &mut solutions {
	// 	*s = Control::random();
	// }

	let mut best_fitness = 9999999999.;
	loop {
		sort_solutions(&mut solutions, &vessel);

		if solutions[0].fitness() < best_fitness {
			best_fitness = solutions[0].fitness();
			println!("New best fitness: {}", best_fitness);
			println!("{:?}", solutions[0]);

			plot.add_trajectory(&vessel, &solutions[0]);
			plot.refresh();

			// if best_fitness < 100. {
			// 	break;
			// }
		}

		//kill the lowest performing half
		solutions.truncate(size/2);

		//reproduce everyone
		for i in 0..size/2 {
			let spawn = solutions[i].mutated();
			solutions.push(spawn);
		}
	}
}
