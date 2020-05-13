% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh3m2DE2
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
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh3m2DE2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_jacobiR_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2DE2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_jacobiR_rot_sym_varpar: pkin has to be [18x1] (double)');
JR_rot=NaN(9,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:22
	% EndTime: 2020-05-07 02:13:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:23
	% EndTime: 2020-05-07 02:13:23
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
	% StartTime: 2020-05-07 02:13:23
	% EndTime: 2020-05-07 02:13:23
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
	% StartTime: 2020-05-07 02:13:24
	% EndTime: 2020-05-07 02:13:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (20->5), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t34 = cos(qJ(1));
	t33 = sin(qJ(1));
	t32 = qJ(2) + qJ(3);
	t31 = cos(t32);
	t30 = sin(t32);
	t29 = t34 * t31;
	t28 = t34 * t30;
	t27 = t33 * t31;
	t26 = t33 * t30;
	t1 = [t27, t28, t28, 0; -t29, t26, t26, 0; 0, -t31, -t31, 0; -t26, t29, t29, 0; t28, t27, t27, 0; 0, t30, t30, 0; t34, 0, 0, 0; t33, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:24
	% EndTime: 2020-05-07 02:13:25
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (2464->13), mult. (4296->20), div. (72->0), fcn. (6486->15), ass. (0->23)
	t88 = sin(pkin(16));
	t89 = cos(pkin(16));
	t93 = sin(pkin(15));
	t97 = cos(pkin(15));
	t83 = t88 * t97 + t89 * t93;
	t84 = -t88 * t93 + t89 * t97;
	t90 = sin(qJ(3));
	t94 = cos(qJ(3));
	t100 = -t83 * t94 - t90 * t84;
	t91 = sin(qJ(2));
	t95 = cos(qJ(2));
	t98 = t90 * t83 - t84 * t94;
	t111 = -t91 * t100 + t98 * t95;
	t99 = t100 * t95 + t91 * t98;
	t96 = cos(qJ(1));
	t92 = sin(qJ(1));
	t87 = pkin(17) + pkin(18);
	t86 = cos(t87);
	t85 = sin(t87);
	t68 = qJ(2) + qJ(3) + atan2(t111 * t85 + t86 * t99, -t111 * t86 + t99 * t85);
	t67 = cos(t68);
	t66 = sin(t68);
	t1 = [t92 * t67, 0, 0, 0; -t96 * t67, 0, 0, 0; 0, 0, 0, 0; -t92 * t66, 0, 0, 0; t96 * t66, 0, 0, 0; 0, 0, 0, 0; t96, 0, 0, 0; t92, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:31
	% EndTime: 2020-05-07 02:13:32
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (3927->18), mult. (6838->28), div. (108->0), fcn. (10332->17), ass. (0->33)
	t268 = sin(qJ(2));
	t273 = cos(qJ(2));
	t264 = sin(pkin(16));
	t265 = cos(pkin(16));
	t270 = sin(pkin(15));
	t275 = cos(pkin(15));
	t259 = t264 * t275 + t265 * t270;
	t260 = -t264 * t270 + t265 * t275;
	t267 = sin(qJ(3));
	t272 = cos(qJ(3));
	t276 = t267 * t259 - t260 * t272;
	t278 = -t259 * t272 - t267 * t260;
	t297 = -t268 * t278 + t276 * t273;
	t277 = t268 * t276 + t273 * t278;
	t266 = sin(qJ(4));
	t269 = sin(qJ(1));
	t286 = t269 * t266;
	t271 = cos(qJ(4));
	t285 = t269 * t271;
	t274 = cos(qJ(1));
	t284 = t274 * t266;
	t283 = t274 * t271;
	t263 = pkin(17) + pkin(18);
	t262 = cos(t263);
	t261 = sin(t263);
	t244 = qJ(2) + qJ(3) + atan2(t261 * t297 + t262 * t277, t277 * t261 - t297 * t262);
	t243 = cos(t244);
	t242 = sin(t244);
	t241 = t243 * t283 - t286;
	t240 = t243 * t284 + t285;
	t239 = t243 * t285 + t284;
	t238 = t243 * t286 - t283;
	t1 = [t239, 0, 0, t240; -t241, 0, 0, t238; 0, 0, 0, t242 * t266; -t238, 0, 0, t241; t240, 0, 0, t239; 0, 0, 0, t242 * t271; t269 * t242, 0, 0, 0; -t274 * t242, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:52
	% EndTime: 2020-05-07 02:13:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (34->8), mult. (68->12), div. (0->0), fcn. (110->8), ass. (0->17)
	t28 = sin(pkin(15));
	t29 = sin(pkin(14));
	t32 = cos(pkin(15));
	t33 = cos(pkin(14));
	t24 = t28 * t33 - t32 * t29;
	t25 = t28 * t29 + t32 * t33;
	t26 = sin(qJ(2));
	t30 = cos(qJ(2));
	t20 = -t30 * t24 - t26 * t25;
	t27 = sin(qJ(1));
	t37 = t20 * t27;
	t21 = -t26 * t24 + t25 * t30;
	t36 = t21 * t27;
	t31 = cos(qJ(1));
	t35 = t31 * t20;
	t34 = t31 * t21;
	t1 = [-t36, t35, 0, 0; t34, t37, 0, 0; 0, t21, 0, 0; -t37, -t34, 0, 0; t35, -t36, 0, 0; 0, t20, 0, 0; t31, 0, 0, 0; t27, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiR_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:52
	% EndTime: 2020-05-07 02:13:52
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (245->15), mult. (458->21), div. (36->2), fcn. (710->11), ass. (0->23)
	t37 = sin(pkin(18));
	t38 = cos(pkin(18));
	t41 = sin(pkin(15));
	t44 = cos(pkin(15));
	t35 = t44 * t37 + t41 * t38;
	t36 = -t41 * t37 + t44 * t38;
	t39 = sin(qJ(2));
	t42 = cos(qJ(2));
	t33 = t39 * t35 - t42 * t36;
	t34 = t42 * t35 + t39 * t36;
	t49 = t34 ^ 2 / t33 ^ 2;
	t29 = qJ(2) + atan2(t34, t33);
	t27 = sin(t29);
	t40 = sin(qJ(1));
	t48 = t40 * t27;
	t28 = cos(t29);
	t47 = t40 * t28;
	t43 = cos(qJ(1));
	t46 = t43 * t27;
	t45 = t43 * t28;
	t30 = 0.1e1 / (0.1e1 + t49);
	t26 = -t30 * t49 - t30 + 0.1e1;
	t1 = [-t47, -t26 * t46, 0, 0; t45, -t26 * t48, 0, 0; 0, t26 * t28, 0, 0; t48, -t26 * t45, 0, 0; -t46, -t26 * t47, 0, 0; 0, -t26 * t27, 0, 0; t43, 0, 0, 0; t40, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiR_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:54
	% EndTime: 2020-05-07 02:13:54
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (1620->25), mult. (2572->33), div. (108->2), fcn. (3942->16), ass. (0->32)
	t103 = sin(qJ(2));
	t107 = cos(qJ(2));
	t105 = sin(pkin(15));
	t106 = cos(qJ(3));
	t109 = cos(pkin(15));
	t120 = sin(qJ(3));
	t95 = t106 * t105 + t120 * t109;
	t96 = -t120 * t105 + t106 * t109;
	t110 = t103 * t96 + t95 * t107;
	t87 = -t103 * t95 + t107 * t96;
	t101 = sin(pkin(18));
	t102 = cos(pkin(18));
	t93 = t109 * t101 + t105 * t102;
	t94 = -t105 * t101 + t109 * t102;
	t85 = t103 * t93 - t107 * t94;
	t86 = t103 * t94 + t107 * t93;
	t118 = t86 ^ 2 / t85 ^ 2;
	t104 = sin(qJ(1));
	t100 = pkin(17) + pkin(18);
	t98 = sin(t100);
	t99 = cos(t100);
	t75 = -qJ(2) - atan2(t86, t85) + pkin(17) - atan2(-t110 * t99 - t98 * t87, t110 * t98 - t87 * t99);
	t73 = sin(t75);
	t114 = t104 * t73;
	t74 = cos(t75);
	t113 = t104 * t74;
	t108 = cos(qJ(1));
	t112 = t108 * t73;
	t111 = t108 * t74;
	t81 = 0.1e1 / (0.1e1 + t118);
	t71 = t81 * t118 + t81 - 0.2e1;
	t1 = [t113, t71 * t112, -t112, 0; -t111, t71 * t114, -t114, 0; 0, t71 * t74, -t74, 0; t114, -t71 * t111, t111, 0; -t112, -t71 * t113, t113, 0; 0, t71 * t73, -t73, 0; t108, 0, 0, 0; t104, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
end