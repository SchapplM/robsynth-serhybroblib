% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh2m2OL_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m2OL_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_jacobiR_rot_sym_varpar: pkin has to be [5x1] (double)');
JR_rot=NaN(9,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t8 = sin(qJ(2));
	t9 = sin(qJ(1));
	t15 = t9 * t8;
	t11 = cos(qJ(1));
	t14 = t11 * t8;
	t10 = cos(qJ(2));
	t13 = t9 * t10;
	t12 = t11 * t10;
	t1 = [-t13, -t14, 0, 0, 0, 0; t12, -t15, 0, 0, 0, 0; 0, t10, 0, 0, 0, 0; t15, -t12, 0, 0, 0, 0; -t14, -t13, 0, 0, 0, 0; 0, -t8, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t20 = qJ(2) + qJ(3);
	t18 = sin(t20);
	t21 = sin(qJ(1));
	t26 = t21 * t18;
	t19 = cos(t20);
	t25 = t21 * t19;
	t22 = cos(qJ(1));
	t24 = t22 * t18;
	t23 = t22 * t19;
	t1 = [-t25, -t24, -t24, 0, 0, 0; t23, -t26, -t26, 0, 0, 0; 0, t19, t19, 0, 0, 0; t26, -t23, -t23, 0, 0, 0; -t24, -t25, -t25, 0, 0, 0; 0, -t18, -t18, 0, 0, 0; t22, 0, 0, 0, 0, 0; t21, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (61->18), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t26 = qJ(2) + qJ(3) + qJ(4);
	t24 = sin(t26);
	t27 = sin(qJ(1));
	t32 = t27 * t24;
	t25 = cos(t26);
	t31 = t27 * t25;
	t28 = cos(qJ(1));
	t30 = t28 * t24;
	t29 = t28 * t25;
	t1 = [-t31, -t30, -t30, -t30, 0, 0; t29, -t32, -t32, -t32, 0, 0; 0, t25, t25, t25, 0, 0; t32, -t29, -t29, -t29, 0, 0; -t30, -t31, -t31, -t31, 0, 0; 0, -t24, -t24, -t24, 0, 0; t28, 0, 0, 0, 0, 0; t27, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (106->23), mult. (20->4), div. (0->0), fcn. (50->4), ass. (0->10)
	t32 = qJ(2) + qJ(3) + qJ(4) + qJ(5);
	t30 = sin(t32);
	t33 = sin(qJ(1));
	t38 = t33 * t30;
	t31 = cos(t32);
	t37 = t33 * t31;
	t34 = cos(qJ(1));
	t36 = t34 * t30;
	t35 = t34 * t31;
	t1 = [-t37, -t36, -t36, -t36, -t36, 0; t35, -t38, -t38, -t38, -t38, 0; 0, t31, t31, t31, t31, 0; t38, -t35, -t35, -t35, -t35, 0; -t36, -t37, -t37, -t37, -t37, 0; 0, -t30, -t30, -t30, -t30, 0; t34, 0, 0, 0, 0, 0; t33, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:58
	% EndTime: 2020-05-03 06:34:58
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (179->34), mult. (76->20), div. (0->0), fcn. (132->6), ass. (0->24)
	t95 = qJ(2) + qJ(3) + qJ(4) + qJ(5);
	t94 = cos(t95);
	t96 = sin(qJ(6));
	t108 = t94 * t96;
	t97 = sin(qJ(1));
	t107 = t97 * t94;
	t106 = t97 * t96;
	t98 = cos(qJ(6));
	t105 = t97 * t98;
	t99 = cos(qJ(1));
	t104 = t99 * t94;
	t103 = t99 * t96;
	t102 = t99 * t98;
	t93 = sin(t95);
	t101 = t93 * t105;
	t100 = t93 * t102;
	t92 = t94 * t98;
	t91 = t93 * t103;
	t90 = t93 * t106;
	t89 = t94 * t102 - t106;
	t88 = -t94 * t103 - t105;
	t87 = -t94 * t105 - t103;
	t86 = t94 * t106 - t102;
	t1 = [t87, -t100, -t100, -t100, -t100, t88; t89, -t101, -t101, -t101, -t101, -t86; 0, t92, t92, t92, t92, -t93 * t96; t86, t91, t91, t91, t91, -t89; t88, t90, t90, t90, t90, t87; 0, -t108, -t108, -t108, -t108, -t93 * t98; t97 * t93, -t104, -t104, -t104, -t104, 0; -t99 * t93, -t107, -t107, -t107, -t107, 0; 0, -t93, -t93, -t93, -t93, 0;];
	JR_rot = t1;
end