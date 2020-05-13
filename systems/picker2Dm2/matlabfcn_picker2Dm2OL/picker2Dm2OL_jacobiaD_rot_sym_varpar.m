% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% picker2Dm2OL
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in picker2Dm2OL_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% 
% Output:
% JaD_rot [3x12]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = picker2Dm2OL_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_jacobiaD_rot_sym_varpar: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm2OL_jacobiaD_rot_sym_varpar: qJD has to be [12x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'picker2Dm2OL_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
JaD_rot=NaN(3,12);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (74->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (195->0), mult. (111->0), div. (45->0), fcn. (66->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (195->0), mult. (111->0), div. (45->0), fcn. (66->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (195->0), mult. (111->0), div. (45->0), fcn. (66->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (74->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiaD_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (372->0), mult. (148->0), div. (60->0), fcn. (88->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiaD_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (372->0), mult. (148->0), div. (60->0), fcn. (88->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
end