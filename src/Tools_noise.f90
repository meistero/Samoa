! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

MODULE Tools_noise
	use SFC_data_types

    implicit none

	integer (kind = selected_int_kind(4)), dimension(0:511), parameter :: p = [ &
	151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,36,103, &
	30,69,142,8,99,37,240,21,10,23,190,6,148,247,120,234,75,0,26,197, &
	62,94,252,219,203,117,35,11,32,57,177,33,88,237,149,56,87,174,20, &
	125,136,171,168,68,175,74,165,71,134,139,48,27,166,77,146,158,231, &
	83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,102, &
	143,54,65,25,63,161,1,216,80,73,209,76,132,187,208,89,18,169,200, &
	196,135,130,116,188,159,86,164,100,109,198,173,186,3,64,52,217,226, &
	250,124,123,5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16, &
	58,17,182,189,28,42,223,183,170,213,119,248,152,2,44,154,163,70, &
	221,153,101,155,167,43,172,9,129,22,39,253,19,98,108,110,79,113, &
	224,232,178,185,112,104,218,246,97,228,251,34,242,193,238,210,144, &
	12,191,179,162,241,81,51,145,235,249,14,239,107,49,192,214,31,181, &
	199,106,157,184,84,204,176,115,121,50,45,127,4,150,254,138,236, &
	205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180, &
	151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,36,103, &
	30,69,142,8,99,37,240,21,10,23,190,6,148,247,120,234,75,0,26,197, &
	62,94,252,219,203,117,35,11,32,57,177,33,88,237,149,56,87,174,20, &
	125,136,171,168,68,175,74,165,71,134,139,48,27,166,77,146,158,231, &
	83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,102, &
	143,54,65,25,63,161,1,216,80,73,209,76,132,187,208,89,18,169,200, &
	196,135,130,116,188,159,86,164,100,109,198,173,186,3,64,52,217,226, &
	250,124,123,5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16, &
	58,17,182,189,28,42,223,183,170,213,119,248,152,2,44,154,163,70, &
	221,153,101,155,167,43,172,9,129,22,39,253,19,98,108,110,79,113, &
	224,232,178,185,112,104,218,246,97,228,251,34,242,193,238,210,144, &
	12,191,179,162,241,81,51,145,235,249,14,239,107,49,192,214,31,181, &
	199,106,157,184,84,204,176,115,121,50,45,127,4,150,254,138,236, &
	205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180 ]

	PRIVATE
	PUBLIC t_noise_1D, t_noise_2D, t_noise_3D

	contains

	pure function t_noise_1D(x, i_depth, r_roughness) result(r_noise)
		real (kind = GRID_SR), intent(in)					:: x				!< position in world coordinates
		integer (kind = GRID_SI), intent(in)				:: i_depth			!< number of recursive noise levels to add (irregularity)
		real (kind = GRID_SR), intent(in)					:: r_roughness		!< roughness of the island boundaries, 0.0 = totally smooth, 1.0 = white noise [0.0, 1.0]
		real (kind = GRID_SR)								:: r_noise			!< resulting noise value

		real (kind = GRID_SR)								:: r_x
		real (kind = GRID_SR)								:: r_factor
		real (kind = GRID_SR)								:: r_sum_factor
		integer (kind = GRID_SI)							:: i_level

		r_noise = 0.0_GRID_SR
		r_factor = 1.0_GRID_SR
		r_sum_factor = 0.0_GRID_SR
		r_x = x

		do i_level = 0, i_depth
			r_noise = r_noise + r_factor * pnoise_1D(r_x)
			r_sum_factor = r_sum_factor + r_factor
			r_factor = r_factor * r_roughness
			r_x = r_x * 2.0_GRID_SR
		end do

		r_noise = r_noise / r_sum_factor
	end function

	pure function t_noise_2D(x, i_depth, r_roughness) result(r_noise)
		real (kind = GRID_SR), dimension(:), intent(in)		:: x				!< position in world coordinates
		integer (kind = GRID_SI), intent(in)				:: i_depth			!< number of recursive noise levels to add (irregularity)
		real (kind = GRID_SR), intent(in)					:: r_roughness		!< roughness of the island boundaries, 0.0 = totally smooth, 1.0 = white noise [0.0, 1.0]
		real (kind = GRID_SR)								:: r_noise			!< resulting noise value

		real (kind = GRID_SR), dimension(2)					:: r_x
		real (kind = GRID_SR)								:: r_factor
		real (kind = GRID_SR)								:: r_sum_factor
		integer (kind = GRID_SI)							:: i_level

		r_noise = 0.0_GRID_SR
		r_factor = 1.0_GRID_SR
		r_sum_factor = 0.0_GRID_SR
		r_x = x

		do i_level = 0, i_depth
			r_noise = r_noise + r_factor * pnoise_2D(r_x)
			r_sum_factor = r_sum_factor + r_factor
			r_factor = r_factor * r_roughness
			r_x = r_x * 2.0_GRID_SR
		end do

		r_noise = r_noise / r_sum_factor
	end function

	pure function t_noise_3D(x, i_depth, r_roughness) result(r_noise)
		real (kind = GRID_SR), dimension(:), intent(in)		:: x				!< position in world coordinates
		integer (kind = GRID_SI), intent(in)				:: i_depth			!< number of recursive noise levels to add (irregularity)
		real (kind = GRID_SR), intent(in)					:: r_roughness		!< roughness of the island boundaries, 0.0 = totally smooth, 1.0 = white noise [0.0, 1.0]
		real (kind = GRID_SR)								:: r_noise			!< resulting noise value

		real (kind = GRID_SR), dimension(3)					:: r_x
		real (kind = GRID_SR)								:: r_factor
		real (kind = GRID_SR)								:: r_sum_factor
		integer (kind = GRID_SI)							:: i_level

		r_noise = 0.0_GRID_SR
		r_factor = 1.0_GRID_SR
		r_sum_factor = 0.0_GRID_SR
		r_x = x

		do i_level = 0, i_depth
			r_noise = r_noise + r_factor * pnoise_3D(r_x)
			r_sum_factor = r_sum_factor + r_factor
			r_factor = r_factor * r_roughness
			r_x = r_x * 2.0_GRID_SR
		end do

		r_noise = r_noise / r_sum_factor
	end function

	pure function lerp(t, a, b)
		real (kind = GRID_SR), intent(in)			:: t			!< interpolator
		real (kind = GRID_SR), intent(in)			:: a, b			!< position in world coordinates
		real (kind = GRID_SR)						:: lerp

    	lerp = a + t * (b - a)
	end function

	elemental function fade(t)
		real (kind = GRID_SR), intent(in)			:: t
		real (kind = GRID_SR)						:: fade

	    fade = t * t * t * (t * (t * 6.0_GRID_SR - 15.0_GRID_SR) + 10.0_GRID_SR)
	end function

	pure function grad_1D(hash, x) result(grad)
	    integer (kind = GRID_SI), intent(in)				:: hash
	    real (kind = GRID_SR), intent(in)					:: x
		real (kind = GRID_SR)								:: grad

	    integer (kind = GRID_SI)							:: h
	    real (kind = GRID_SR), dimension(2)					:: u

	    h = iand(hash, 15)

	    if (h < 8) then
	        u(1) = x
	    else
	        u(1) = 0
	    end if

	    if (h < 4) then
	        u(2) = 0
	    elseif (h == 12 .or. h == 14) then
	        u(2) = x
	    else
	        u(2) = 0
	    end if

	    if (iand(h, 1) .ne. 0) then
	        u(1) = -u(1)
	    end if

	    if (iand(h, 2) .ne. 0) then
	        u(2) = -u(2)
	    end if

	    grad = u(1) + u(2)
	end function

	pure function grad_2D(hash, x) result(grad)
	    integer (kind = GRID_SI), intent(in)				:: hash
	    real (kind = GRID_SR), dimension(:), intent(in)		:: x
		real (kind = GRID_SR)								:: grad

	    integer (kind = GRID_SI)							:: h
	    real (kind = GRID_SR), dimension(2)					:: u

	    h = iand(hash, 15)

	    if (h < 8) then
	        u(1) = x(1)
	    else
	        u(1) = x(2)
	    end if

	    if (h < 4) then
	        u(2) = x(2)
	    elseif (h == 12 .or. h == 14) then
	        u(2) = x(1)
	    else
	        u(2) = 0
	    end if

	    if (iand(h, 1) .ne. 0) then
	        u(1) = -u(1)
	    end if

	    if (iand(h, 2) .ne. 0) then
	        u(2) = -u(2)
	    end if

	    grad = u(1) + u(2)
	end function

	pure function grad_3D(hash, x) result(grad)
	    integer (kind = GRID_SI), intent(in)				:: hash
	    real (kind = GRID_SR), dimension(:), intent(in)		:: x
		real (kind = GRID_SR)								:: grad

	    integer (kind = GRID_SI)							:: h
	    real (kind = GRID_SR), dimension(2)					:: u

	    h = iand(hash, 15)

	    if (h < 8) then
	        u(1) = x(1)
	    else
	        u(1) = x(2)
	    end if

	    if (h < 4) then
	        u(2) = x(2)
	    elseif (h == 12 .or. h == 14) then
	        u(2) = x(1)
	    else
	        u(2) = x(3)
	    end if

	    if (iand(h, 1) .ne. 0) then
	        u(1) = -u(1)
	    end if

	    if (iand(h, 2) .ne. 0) then
	        u(2) = -u(2)
	    end if

	    grad = u(1) + u(2)
	end function

	pure function pnoise_1D(x) result(pnoise)
	    real (kind = GRID_SR), intent(in)				:: x
	    real (kind = GRID_SR)							:: pnoise

	    real (kind = GRID_SR) 							:: r_x
	    real (kind = GRID_SR) 							:: u
	    integer (kind = GRID_SI)						:: i_x
	    integer (kind = GRID_SI)						:: i_hash_A, i_hash_B

	    real (kind = GRID_SR)							:: r_grad_A, r_grad_B


	    i_x = iand(int(floor(x)), 255)

	    i_hash_A =  p(i_x)
	    i_hash_B =  p(iand((i_x + 1), 255))

	    r_x = x - floor(x)

	    r_grad_A =  grad_1D(i_hash_A, r_x)
	    r_grad_B =  grad_1D(i_hash_B, r_x - 1.0_GRID_SR)

	    u = fade(r_x)
	    pnoise = lerp(u, r_grad_A, r_grad_B)
	end function

	pure function pnoise_2D(x) result(pnoise)
	    real (kind = GRID_SR), dimension(:), intent(in)	:: x
	    real (kind = GRID_SR)							:: pnoise

	    real (kind = GRID_SR), dimension(2)				:: r_x
	    real (kind = GRID_SR), dimension(2)				:: u
	    integer (kind = GRID_SI), dimension(2)			:: i_x
	    integer (kind = GRID_SI)						:: i_hash_A, i_hash_B, i_hash_AA, i_hash_AB, i_hash_BA, i_hash_BB

	    real (kind = GRID_SR)							:: r_grad_AA, r_grad_AB, r_grad_BA, r_grad_BB

	    i_x = iand(int(floor(x)), 255)

	    i_hash_A =  p(i_x(1))
	    i_hash_B =  p(iand((i_x(1) + 1), 255))

	    i_hash_AA = p(i_hash_A + i_x(2))
	    i_hash_BA = p(i_hash_B + i_x(2))
	    i_hash_AB = p(iand((i_hash_A + i_x(2) + 1), 255))
	    i_hash_BB = p(iand((i_hash_B + i_x(2) + 1), 255))

	    r_x = x - floor(x)

	    r_grad_AA =  grad_2D(i_hash_AA, r_x)
	    r_grad_BA =  grad_2D(i_hash_BA, r_x - [ 1.0_GRID_SR, 0.0_GRID_SR ])
	    r_grad_AB =  grad_2D(i_hash_AB, r_x - [ 0.0_GRID_SR, 1.0_GRID_SR ])
	    r_grad_BB =  grad_2D(i_hash_BB, r_x - [ 1.0_GRID_SR, 1.0_GRID_SR ])

	    u = fade(r_x)
	    pnoise = lerp(u(2), lerp(u(1), r_grad_AA, r_grad_BA), lerp(u(1), r_grad_AB, r_grad_BB))
	end function

	pure function pnoise_3D(x) result(pnoise)
	    real (kind = GRID_SR), dimension(:), intent(in)	:: x
	    real (kind = GRID_SR)							:: pnoise

	    real (kind = GRID_SR), dimension(3)				:: r_x
	    real (kind = GRID_SR), dimension(3)				:: u
	    integer (kind = GRID_SI), dimension(3)			:: i_x
	    integer (kind = GRID_SI)						:: i_hash_A, i_hash_B, i_hash_AA, i_hash_AB, i_hash_BA, i_hash_BB, &
	    													i_hash_AAA, i_hash_AAB, i_hash_ABA, i_hash_ABB, &
	    													i_hash_BAA, i_hash_BAB, i_hash_BBA, i_hash_BBB

	    real (kind = GRID_SR)							:: r_grad_AAA, r_grad_AAB, r_grad_ABA, r_grad_ABB, &
	    													r_grad_BAA, r_grad_BAB, r_grad_BBA, r_grad_BBB


	    i_x = iand(int(floor(x)), 255)
	    i_hash_A =  p(i_x(1))
	    i_hash_B =  p(iand((i_x(1) + 1), 255))

	    i_hash_AA = p(i_hash_A + i_x(2))
	    i_hash_BA = p(i_hash_B + i_x(2))
	    i_hash_AB = p(iand((i_hash_A + i_x(2) + 1), 255))
	    i_hash_BB = p(iand((i_hash_B + i_x(2) + 1), 255))

	    i_hash_AAA = p(i_hash_AA + i_x(3))
	    i_hash_ABA = p(i_hash_AB + i_x(3))
	    i_hash_BAA = p(i_hash_BA + i_x(3))
	    i_hash_BBA = p(i_hash_BB + i_x(3))
	    i_hash_AAB = p(i_hash_AA + i_x(3) + 1)
	    i_hash_BAB = p(i_hash_BA + i_x(3) + 1)
	    i_hash_ABB = p(i_hash_AB + i_x(3) + 1)
	    i_hash_BBB = p(i_hash_BB + i_x(3) + 1)

	    r_x = x - floor(x)

	    r_grad_AAA = grad_3D(i_hash_AAA, r_x)
	    r_grad_BAA = grad_3D(i_hash_BAA, r_x - [ 1.0_GRID_SR, 0.0_GRID_SR, 0.0_GRID_SR ])
	    r_grad_ABA = grad_3D(i_hash_ABA, r_x - [ 0.0_GRID_SR, 1.0_GRID_SR, 0.0_GRID_SR ])
	    r_grad_BBA = grad_3D(i_hash_BBA, r_x - [ 1.0_GRID_SR, 1.0_GRID_SR, 0.0_GRID_SR ])
	    r_grad_AAB = grad_3D(i_hash_AAB, r_x - [ 0.0_GRID_SR, 0.0_GRID_SR, 1.0_GRID_SR ])
	    r_grad_BAB = grad_3D(i_hash_BAB, r_x - [ 1.0_GRID_SR, 0.0_GRID_SR, 1.0_GRID_SR ])
	    r_grad_ABB = grad_3D(i_hash_ABB, r_x - [ 0.0_GRID_SR, 1.0_GRID_SR, 1.0_GRID_SR ])
	    r_grad_BBB = grad_3D(i_hash_BBB, r_x - [ 1.0_GRID_SR, 1.0_GRID_SR, 1.0_GRID_SR ])

	    u = fade(r_x)
	    pnoise = lerp(u(3), lerp(u(2), lerp(u(1), r_grad_AAA, r_grad_BAA), lerp(u(1), r_grad_ABA, r_grad_BBA)), lerp(u(2), lerp(u(1), r_grad_AAB, r_grad_BAB), lerp(u(1), r_grad_ABB, r_grad_BBB)))
	end function
end MODULE
