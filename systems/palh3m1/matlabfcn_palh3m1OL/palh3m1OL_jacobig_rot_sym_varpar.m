% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh3m1OL
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% Jg_rot [3x10]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:16
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = palh3m1OL_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),uint8(0),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1OL_jacobig_rot_sym_varpar: qJ has to be [10x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1OL_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1OL_jacobig_rot_sym_varpar: pkin has to be [16x1] (double)');
Jg_rot=NaN(3,10);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:35
	% EndTime: 2020-04-20 17:16:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:35
	% EndTime: 2020-04-20 17:16:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:35
	% EndTime: 2020-04-20 17:16:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)), 0, 0, 0, 0, 0, 0, 0, 0; 0, -cos(qJ(1)), 0, 0, 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:35
	% EndTime: 2020-04-20 17:16:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t25 = cos(qJ(1));
	t24 = sin(qJ(1));
	t1 = [0, t24, t24, 0, 0, 0, 0, 0, 0, 0; 0, -t25, -t25, 0, 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:35
	% EndTime: 2020-04-20 17:16:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (6->2), ass. (0->3)
	t30 = cos(qJ(1));
	t29 = sin(qJ(1));
	t1 = [0, t29, t29, t29, 0, 0, 0, 0, 0, 0; 0, -t30, -t30, -t30, 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:36
	% EndTime: 2020-04-20 17:16:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->6), mult. (2->2), div. (0->0), fcn. (11->4), ass. (0->5)
	t88 = cos(qJ(1));
	t87 = sin(qJ(1));
	t86 = qJ(2) + qJ(3) + qJ(4);
	t85 = sin(t86);
	t1 = [0, t87, t87, t87, -t88 * t85, 0, 0, 0, 0, 0; 0, -t88, -t88, -t88, -t87 * t85, 0, 0, 0, 0, 0; 1, 0, 0, 0, cos(t86), 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:35
	% EndTime: 2020-04-20 17:16:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, sin(qJ(1)), 0, 0, 0, 0; 0, 0, 0, 0, 0, -cos(qJ(1)), 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobig_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:35
	% EndTime: 2020-04-20 17:16:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t17 = cos(qJ(1));
	t16 = sin(qJ(1));
	t1 = [0, t16, 0, 0, 0, 0, t16, 0, 0, 0; 0, -t17, 0, 0, 0, 0, -t17, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobig_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 17:16:35
	% EndTime: 2020-04-20 17:16:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (6->2), ass. (0->3)
	t24 = cos(qJ(1));
	t23 = sin(qJ(1));
	t1 = [0, t23, 0, 0, 0, 0, t23, t23, 0, 0; 0, -t24, 0, 0, 0, 0, -t24, -t24, 0, 0; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
end