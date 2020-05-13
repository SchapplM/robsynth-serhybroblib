% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh3m2TE
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% JR_rot [9x4]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh3m2TE_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_jacobiR_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2TE_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_jacobiR_rot_sym_varpar: pkin has to be [18x1] (double)');
JR_rot=NaN(9,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:40
	% EndTime: 2020-05-07 01:41:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:40
	% EndTime: 2020-05-07 01:41:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0; t9, 0, 0, 0; 0, 0, 0, 0; -t9, 0, 0, 0; -t8, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:41
	% EndTime: 2020-05-07 01:41:41
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
	t1 = [-t13, -t14, 0, 0; t12, -t15, 0, 0; 0, t10, 0, 0; t15, -t12, 0, 0; -t14, -t13, 0, 0; 0, -t8, 0, 0; t11, 0, 0, 0; t9, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:41
	% EndTime: 2020-05-07 01:41:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->8), mult. (44->8), div. (0->0), fcn. (78->6), ass. (0->13)
	t41 = cos(qJ(2));
	t40 = cos(qJ(3));
	t36 = sin(qJ(3));
	t37 = sin(qJ(2));
	t33 = t36 * t37 - t40 * t41;
	t38 = sin(qJ(1));
	t29 = t38 * t33;
	t39 = cos(qJ(1));
	t31 = t39 * t33;
	t34 = t36 * t41 + t40 * t37;
	t32 = t39 * t34;
	t30 = t38 * t34;
	t1 = [-t29, t32, t32, 0; t31, t30, t30, 0; 0, t33, t33, 0; -t30, -t31, -t31, 0; t32, -t29, -t29, 0; 0, t34, t34, 0; t39, 0, 0, 0; t38, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:42
	% EndTime: 2020-05-07 01:41:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->7), mult. (28->12), div. (0->0), fcn. (46->8), ass. (0->14)
	t41 = cos(pkin(15));
	t40 = cos(qJ(1));
	t39 = sin(pkin(15));
	t38 = sin(qJ(1));
	t37 = cos(pkin(16));
	t36 = sin(pkin(16));
	t35 = pkin(17) + pkin(18);
	t34 = cos(t35);
	t33 = sin(t35);
	t32 = -t36 * t39 + t37 * t41;
	t31 = t36 * t41 + t37 * t39;
	t30 = t33 * t31 - t32 * t34;
	t29 = t31 * t34 + t32 * t33;
	t1 = [-t30 * t38, 0, 0, 0; t30 * t40, 0, 0, 0; 0, 0, 0, 0; t38 * t29, 0, 0, 0; -t40 * t29, 0, 0, 0; 0, 0, 0, 0; t40, 0, 0, 0; t38, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:43
	% EndTime: 2020-05-07 01:41:43
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (71->14), mult. (116->22), div. (0->0), fcn. (172->10), ass. (0->22)
	t101 = pkin(17) + pkin(18);
	t100 = cos(t101);
	t102 = sin(pkin(16));
	t103 = cos(pkin(16));
	t106 = sin(pkin(15));
	t109 = cos(pkin(15));
	t97 = t102 * t109 + t103 * t106;
	t98 = -t102 * t106 + t103 * t109;
	t99 = sin(t101);
	t112 = t100 * t98 - t97 * t99;
	t105 = sin(qJ(1));
	t111 = t112 * t105;
	t108 = cos(qJ(1));
	t110 = t112 * t108;
	t107 = cos(qJ(4));
	t104 = sin(qJ(4));
	t96 = t97 * t100 + t98 * t99;
	t95 = -t105 * t104 + t107 * t110;
	t94 = t104 * t110 + t105 * t107;
	t93 = t108 * t104 + t107 * t111;
	t92 = t104 * t111 - t108 * t107;
	t1 = [t93, 0, 0, t94; -t95, 0, 0, t92; 0, 0, 0, -t96 * t104; -t92, 0, 0, t95; t94, 0, 0, t93; 0, 0, 0, -t96 * t107; -t105 * t96, 0, 0, 0; t108 * t96, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:44
	% EndTime: 2020-05-07 01:41:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (32->7), mult. (68->12), div. (0->0), fcn. (110->8), ass. (0->17)
	t34 = sin(pkin(15));
	t37 = cos(pkin(14));
	t40 = sin(pkin(14));
	t41 = cos(pkin(15));
	t29 = t37 * t34 - t41 * t40;
	t30 = t40 * t34 + t37 * t41;
	t32 = sin(qJ(2));
	t35 = cos(qJ(2));
	t27 = t32 * t29 - t30 * t35;
	t33 = sin(qJ(1));
	t45 = t33 * t27;
	t36 = cos(qJ(1));
	t44 = t36 * t27;
	t38 = -t35 * t29 - t32 * t30;
	t43 = t33 * t38;
	t42 = t36 * t38;
	t1 = [t45, t42, 0, 0; -t44, t43, 0, 0; 0, -t27, 0, 0; -t43, t44, 0, 0; t42, t45, 0, 0; 0, t38, 0, 0; t36, 0, 0, 0; t33, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiR_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:45
	% EndTime: 2020-05-07 01:41:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->4), mult. (12->8), div. (0->0), fcn. (22->6), ass. (0->9)
	t21 = cos(pkin(15));
	t20 = cos(qJ(1));
	t19 = sin(pkin(15));
	t18 = sin(qJ(1));
	t17 = cos(pkin(18));
	t16 = sin(pkin(18));
	t15 = -t16 * t19 + t17 * t21;
	t14 = t16 * t21 + t17 * t19;
	t1 = [t18 * t15, 0, 0, 0; -t20 * t15, 0, 0, 0; 0, 0, 0, 0; t18 * t14, 0, 0, 0; -t20 * t14, 0, 0, 0; 0, 0, 0, 0; t20, 0, 0, 0; t18, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiR_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:45
	% EndTime: 2020-05-07 01:41:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (20->5), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t39 = cos(qJ(1));
	t38 = sin(qJ(1));
	t37 = qJ(3) + qJ(2);
	t36 = cos(t37);
	t35 = sin(t37);
	t34 = t39 * t36;
	t33 = t39 * t35;
	t32 = t38 * t36;
	t31 = t38 * t35;
	t1 = [t32, t33, t33, 0; -t34, t31, t31, 0; 0, -t36, -t36, 0; -t31, t34, t34, 0; t33, t32, t32, 0; 0, t35, t35, 0; t39, 0, 0, 0; t38, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
end