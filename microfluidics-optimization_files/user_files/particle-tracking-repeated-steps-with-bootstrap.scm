;; REPEATED PARTICLE TRACKING
;; Inertial microfluidic devices often have several repeated segments
;; to modify fluid flow and particle positions within a channel.
;; Meshing and simulating the entire device would be computationally expensive
;; and applying periodic boundary conditions assumes the device is infinitely long.
;; 
;; To simulate particle flow in an inertial microfluidic device with finite repeated steps,
;; we first release particles at the inlet from an initial injection DPM definition.
;; Particles are then sampled at the outlet plane and re-injected in subsequent steps.

;; LOAD THE CASE FILE (otherwise this is strangely unreliable)
;(define root "C:\\Users\\rylab\\Desktop\\microfluidics-opt\\microfluidics-optimization_files")
;(ti-menu-load-string (format #f "!copy \"~a\\user_files\\fluent-case-backup\\FFF-Setup-Output.cas.h5\" \"~a\\dp0\\FFF\\Fluent\\FFF-Setup-Output.cas.h5\"" root root))


;; DEFINITIONS
(define particle_diameter 10e-6)
(define particle_radius (/ particle_diameter 2))
(define root "C:\\Users\\rylab\\Desktop\\microfluidics-opt\\microfluidics-optimization_files")

;; FUNCTION DEFINITIONS
(define (append-elt lst x)
  (append lst (list x)))

(define (average li)
  (/ (apply + li) (length li)))

(define (summation li avg)
  (if (null? li)
    0
    (+  
       (* (- (car li) avg) 
          (- (car li) avg))
    (summation (cdr li) avg))))

(define (sd li)
    (sqrt
        (/  
          (summation li (average li))
          (- (length li) 1))))

(define (process-step-data dpm-filename data-filename) (with-input-from-file dpm-filename 
	(lambda ()
		(define next-line (read))
		(define ypos-array '())
		(define zpos-array '())
		(define xvelocity-array '())

		; loop over the outlet file and read the position of each particle at the exit
		; save the relevent data to an array
		(do ((x 0 (+ x 1))) ((eof-object? next-line))
			(if (> x 1) (begin
				; TODO: make sure that I remove any particles that have collided with a wall
				(set! ypos-array (append-elt ypos-array (list-ref (car next-line) 1)))
				(set! zpos-array (append-elt zpos-array (list-ref (car next-line) 2)))
				(set! xvelocity-array (append-elt xvelocity-array (list-ref (car next-line) 3)))
			))
			(set! next-line (read))
		)

		; write the file for each step
		(define port (open-output-file data-filename #t))
		(write-string (format #f "~a,~a,~a,~a,~a,~a"
			; compute the values of x, y, and velocity
			(average ypos-array) (sd ypos-array)
			(average zpos-array) (sd zpos-array)
			(average xvelocity-array) (sd xvelocity-array)
		) port)
		(newline port)
		(close-output-port port)

		; display the results that will have been written to the file for feedback
		(display (format #f "YPos: ~a +/- ~a\n" (average ypos-array) (sd ypos-array)))
		(display (format #f "ZPos: ~a +/- ~a\n" (average zpos-array) (sd zpos-array)))
		(display (format #f "XVelocity: ~a +/- ~a\n" (average xvelocity-array) (sd xvelocity-array)))
	)
))

; set up the variables
(rp-var-define 'volumetric_flow_rate_ul_per_min_scheme 0.0 'real #f)
(rp-var-define 'channel_width_scheme 0.0 'real #f)
(rp-var-define 'channel_height_scheme 0.0 'real #f)
(rp-var-define 'channel_length_scheme 0.0 'real #f)
(rp-var-define 'notch_height_scheme 0.0 'real #f)
(rp-var-define 'notch_length_scheme 0.0 'real #f)
(rp-var-define 're_scheme 0.0 'real #f)
(rp-var-define 'max_number_of_fluent_iterations_scheme 0.0 'real #f)
(rp-var-define 'number_of_steps_scheme 0.0 'real #f)

; execute the port
(ti-menu-load-string (format #f "/define/user-defined/execute-on-demand \"port_input_vars_to_scheme::libudf\""))

; define a string for saving the files for later processing
(define cw (rpgetvar 'channel_width_scheme))
(define ch (rpgetvar 'channel_height_scheme))
(define cl (rpgetvar 'channel_length_scheme))
(define nh (rpgetvar 'notch_height_scheme))
(define nl (rpgetvar 'notch_length_scheme))
(define re (rpgetvar 're_scheme))
(define ninjections (rpgetvar 'number_of_steps_scheme))
(define niterations (rpgetvar 'max_number_of_fluent_iterations_scheme))
(define q (rpgetvar 'volumetric_flow_rate_ul_per_min_scheme))
(define param_string (format #f "cw~a_ch~a_cl~a_nh~a_nl~a_q~a_re~a" cw ch cl nh nl q re))

; edit the analysis surfaces
; inlet surface
(ti-menu-load-string (format #f "/surface/edit-surface particle-inlet particle-inlet three-points 0 ~a ~a 0 ~a ~a 0 ~a ~a yes no" 
	(- (/ (* 1e-6 ch) 2) particle_radius)
	(- (/ (* 1e-6 cw) 2) particle_radius)
	(* (- (/ (* 1e-6 ch) 2) particle_radius) -1)
	(- (/ (* 1e-6 cw) 2) particle_radius)
	(* (- (/ (* 1e-6 ch) 2) particle_radius) -1)
	(* (- (/ (* 1e-6 cw) 2) particle_radius) -1) ; -1 for full width, 0 for symmetry along xy plane
))
; inlet horizontal line
(ti-menu-load-string (format #f "/surface/edit-surface line-inlet-horizontal line-inlet-horizontal 0 0 ~a 0 0 ~a"
	(* (- (/ (* 1e-6 cw) 2) particle_radius) -1)
	(- (/ (* 1e-6 cw) 2) particle_radius)
))
; outlet horizontal line
(ti-menu-load-string (format #f "/surface/edit-surface line-outlet-horizontal line-outlet-horizontal ~a 0 ~a ~a 0 ~a"
	(* 1e-6 cl)
	(* (- (/ (* 1e-6 cw) 2) particle_radius) -1)
	(* 1e-6 cl)
	(- (/ (* 1e-6 cw) 2) particle_radius)
))
; reinjection reference frame
(ti-menu-load-string (format #f "/define/reference-frames/edit/reinjection-reference-frame origin ~a 0 0 ()" (* -1e-6 cl)))



; ==========================================================
; INITIALIZE
; ==========================================================

; run the initizliation once so that we have some consistency with the data
(ti-menu-load-string (format #f "/solve/initialize/hyb-initialization"))

; run the initialization again so that each subsequent time, we can append "o" to overwrite the previous initialization
(ti-menu-load-string (format #f "/solve/initialize/hyb-initialization o"))
(ti-menu-load-string (format #f "/solve/iterate ~a" niterations))

; copy the output parameters results to the data folder
(ti-menu-load-string (format #f "!del /F /Q ~a\\user_files\\data\\*_~a.*" root param_string))
(ti-menu-load-string (format #f "/define/parameters/output-parameters/write-all-to-file \"~a\\user_files\\data\\sim_params_~a.txt\"" root param_string))
(ti-menu-load-string (format #f "/file/write-profile \"~a\\user_files\\data\\inlet_velocity_profile_~a.prof\" line-inlet-horizontal () velocity-magnitude ()" root param_string))
(ti-menu-load-string (format #f "/file/write-profile \"~a\\user_files\\data\\outlet_velocity_profile_~a.prof\" line-outlet-horizontal () velocity-magnitude ()" root param_string))



; ==========================================================
; RUN THE INJECTIONS
; ==========================================================

;; RUN THE INITIAL INJECTION
(display (format #f "\n\n--------------------\nINITIAL INJECTION\n--------------------\n\n"))
(ti-menu-load-string (format #f "/report/dpm-sample initial-injection () outlet () inlet-capture () no"))

; save the relevant output data
(ti-menu-load-string (format #f "!copy \"~a\\dp0\\FFF\\Fluent\\outlet.dpm\" \"~a\\user_files\\data\\dpm\\dpm_data_~a_step~a.dpm\"" root root param_string 1))
(ti-menu-load-string (format #f "!copy \"~a\\dp0\\FFF\\Fluent\\inlet-capture.dpm\" \"~a\\user_files\\data\\dpm\\dpm_inlet_data_~a_step~a.dpm\"" root root param_string 1))
(process-step-data (format #f "~a\\dp0\\FFF\\Fluent\\outlet.dpm" root) (format #f "~a\\user_files\\data\\step_data_~a.csv" root param_string))

;; LOOP THE INJECTIONS
(DO ((x 1 (+ 1 x))) ((> x (- ninjections 1)))
	(ti-menu-load-string (format #f "/define/injections/list-particles file-injection"))

	(display (format #f "\n\n--------------------\nIteration ~a of ~a\n--------------------\n\n" (+ x 1) ninjections))
	(ti-menu-load-string (format #f "/report/dpm-sample file-injection () outlet () inlet-capture () no"))

	; save the relevant output data
	(ti-menu-load-string (format #f "!copy \"~a\\dp0\\FFF\\Fluent\\outlet.dpm\" \"~a\\user_files\\data\\dpm\\dpm_data_~a_step~a.dpm\"" root root param_string (+ x 1)))
	(ti-menu-load-string (format #f "!copy \"~a\\dp0\\FFF\\Fluent\\inlet-capture.dpm\" \"~a\\user_files\\data\\dpm\\dpm_inlet_data_~a_step~a.dpm\"" root root param_string (+ x 1)))
	(process-step-data (format #f "~a\\dp0\\FFF\\Fluent\\outlet.dpm" root) (format #f "~a\\user_files\\data\\step_data_~a.csv" root param_string))
)