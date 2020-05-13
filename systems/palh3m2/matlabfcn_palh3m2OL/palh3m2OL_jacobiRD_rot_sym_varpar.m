% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% JRD_rot [9x10]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = palh3m2OL_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),uint8(0),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_jacobiRD_rot_sym_varpar: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2OL_jacobiRD_rot_sym_varpar: qJD has to be [10x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2OL_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_jacobiRD_rot_sym_varpar: pkin has to be [16x1] (double)');
JRD_rot=NaN(9,10);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t32 = sin(qJ(1));
	t39 = qJD(1) * t32;
	t34 = cos(qJ(1));
	t38 = qJD(1) * t34;
	t31 = sin(qJ(2));
	t37 = qJD(2) * t31;
	t33 = cos(qJ(2));
	t36 = qJD(2) * t33;
	t35 = qJD(2) * t34;
	t30 = t32 * t37 - t33 * t38;
	t29 = t31 * t38 + t32 * t36;
	t28 = t31 * t35 + t33 * t39;
	t27 = t31 * t39 - t33 * t35;
	t1 = [t30, t27, 0, 0, 0, 0, 0, 0, 0, 0; -t28, -t29, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t37, 0, 0, 0, 0, 0, 0, 0, 0; t29, t28, 0, 0, 0, 0, 0, 0, 0, 0; t27, t30, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t36, 0, 0, 0, 0, 0, 0, 0, 0; -t39, 0, 0, 0, 0, 0, 0, 0, 0, 0; t38, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:52
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (57->10), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t77 = qJD(2) + qJD(3);
	t79 = sin(qJ(1));
	t84 = t77 * t79;
	t80 = cos(qJ(1));
	t83 = t77 * t80;
	t82 = qJD(1) * t79;
	t81 = qJD(1) * t80;
	t78 = qJ(2) + qJ(3);
	t76 = cos(t78);
	t75 = sin(t78);
	t74 = t77 * t76;
	t73 = t77 * t75;
	t72 = -t75 * t84 + t76 * t81;
	t71 = t75 * t81 + t76 * t84;
	t70 = t75 * t83 + t76 * t82;
	t69 = -t75 * t82 + t76 * t83;
	t1 = [t72, t69, t69, 0, 0, 0, 0, 0, 0, 0; t70, t71, t71, 0, 0, 0, 0, 0, 0, 0; 0, t73, t73, 0, 0, 0, 0, 0, 0, 0; -t71, -t70, -t70, 0, 0, 0, 0, 0, 0, 0; t69, t72, t72, 0, 0, 0, 0, 0, 0, 0; 0, t74, t74, 0, 0, 0, 0, 0, 0, 0; -t82, 0, 0, 0, 0, 0, 0, 0, 0, 0; t81, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:52
	% EndTime: 2020-05-07 04:44:52
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (137->11), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t94 = qJD(2) + qJD(3) + qJD(4);
	t96 = sin(qJ(1));
	t101 = t94 * t96;
	t97 = cos(qJ(1));
	t100 = t94 * t97;
	t99 = qJD(1) * t96;
	t98 = qJD(1) * t97;
	t95 = qJ(2) + qJ(3) + qJ(4);
	t93 = cos(t95);
	t92 = sin(t95);
	t91 = t94 * t93;
	t90 = t94 * t92;
	t89 = -t92 * t101 + t93 * t98;
	t88 = t93 * t101 + t92 * t98;
	t87 = t92 * t100 + t93 * t99;
	t86 = t93 * t100 - t92 * t99;
	t1 = [t89, t86, t86, t86, 0, 0, 0, 0, 0, 0; t87, t88, t88, t88, 0, 0, 0, 0, 0, 0; 0, t90, t90, t90, 0, 0, 0, 0, 0, 0; -t88, -t87, -t87, -t87, 0, 0, 0, 0, 0, 0; t86, t89, t89, t89, 0, 0, 0, 0, 0, 0; 0, t91, t91, t91, 0, 0, 0, 0, 0, 0; -t99, 0, 0, 0, 0, 0, 0, 0, 0, 0; t98, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:53
	% EndTime: 2020-05-07 04:44:53
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (343->32), mult. (279->57), div. (0->0), fcn. (279->6), ass. (0->41)
	t335 = sin(qJ(1));
	t333 = qJ(2) + qJ(3) + qJ(4);
	t331 = cos(t333);
	t340 = qJD(1) * t331 + qJD(5);
	t330 = sin(t333);
	t332 = qJD(2) + qJD(3) + qJD(4);
	t337 = cos(qJ(1));
	t352 = t332 * t337;
	t343 = t330 * t352;
	t357 = t340 * t335 + t343;
	t356 = t332 * t331;
	t334 = sin(qJ(5));
	t355 = t332 * t334;
	t354 = t332 * t335;
	t336 = cos(qJ(5));
	t353 = t332 * t336;
	t351 = qJD(1) * t335;
	t350 = qJD(1) * t337;
	t349 = qJD(5) * t334;
	t348 = qJD(5) * t336;
	t347 = qJD(5) * t337;
	t346 = t330 * t353;
	t345 = t330 * t354;
	t344 = t331 * t354;
	t342 = t331 * t352;
	t341 = qJD(5) * t331 + qJD(1);
	t339 = t341 * t337;
	t338 = t340 * t337;
	t329 = -t331 * t350 + t345;
	t328 = t331 * t351 + t343;
	t327 = t331 * t349 + t346;
	t326 = -t330 * t355 + t331 * t348;
	t325 = t336 * t344 + (-t335 * t349 + t336 * t350) * t330;
	t324 = -t334 * t344 + (-t334 * t350 - t335 * t348) * t330;
	t323 = t336 * t342 + (-t334 * t347 - t336 * t351) * t330;
	t322 = -t334 * t342 + (t334 * t351 - t336 * t347) * t330;
	t321 = t336 * t338 + (-t341 * t334 - t346) * t335;
	t320 = t341 * t336 * t335 + (t338 - t345) * t334;
	t319 = t334 * t339 + t357 * t336;
	t318 = -t357 * t334 + t336 * t339;
	t1 = [t321, t323, t323, t323, t318, 0, 0, 0, 0, 0; t319, t325, t325, t325, t320, 0, 0, 0, 0, 0; 0, t327, t327, t327, t330 * t348 + t331 * t355, 0, 0, 0, 0, 0; -t320, t322, t322, t322, -t319, 0, 0, 0, 0, 0; t318, t324, t324, t324, t321, 0, 0, 0, 0, 0; 0, t326, t326, t326, -t330 * t349 + t331 * t353, 0, 0, 0, 0, 0; t330 * t350 + t344, t328, t328, t328, 0, 0, 0, 0, 0, 0; t330 * t351 - t342, t329, t329, t329, 0, 0, 0, 0, 0, 0; 0, -t356, -t356, -t356, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t32 = sin(qJ(1));
	t39 = qJD(1) * t32;
	t34 = cos(qJ(1));
	t38 = qJD(1) * t34;
	t31 = sin(qJ(6));
	t37 = qJD(6) * t31;
	t33 = cos(qJ(6));
	t36 = qJD(6) * t33;
	t35 = qJD(6) * t34;
	t30 = t32 * t37 - t33 * t38;
	t29 = t31 * t38 + t32 * t36;
	t28 = t31 * t35 + t33 * t39;
	t27 = t31 * t39 - t33 * t35;
	t1 = [t30, 0, 0, 0, 0, t27, 0, 0, 0, 0; -t28, 0, 0, 0, 0, -t29, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t37, 0, 0, 0, 0; t29, 0, 0, 0, 0, t28, 0, 0, 0, 0; t27, 0, 0, 0, 0, t30, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t36, 0, 0, 0, 0; -t39, 0, 0, 0, 0, 0, 0, 0, 0, 0; t38, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiRD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:52
	% EndTime: 2020-05-07 04:44:52
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (61->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t67 = qJ(2) + qJ(7);
	t64 = sin(t67);
	t66 = qJD(2) + qJD(7);
	t75 = t66 * t64;
	t65 = cos(t67);
	t74 = t66 * t65;
	t68 = sin(qJ(1));
	t73 = t66 * t68;
	t69 = cos(qJ(1));
	t72 = t66 * t69;
	t71 = qJD(1) * t68;
	t70 = qJD(1) * t69;
	t63 = t64 * t73 - t65 * t70;
	t62 = t64 * t70 + t65 * t73;
	t61 = t64 * t72 + t65 * t71;
	t60 = t64 * t71 - t65 * t72;
	t1 = [t63, t60, 0, 0, 0, 0, t60, 0, 0, 0; -t61, -t62, 0, 0, 0, 0, -t62, 0, 0, 0; 0, -t75, 0, 0, 0, 0, -t75, 0, 0, 0; t62, t61, 0, 0, 0, 0, t61, 0, 0, 0; t60, t63, 0, 0, 0, 0, t63, 0, 0, 0; 0, -t74, 0, 0, 0, 0, -t74, 0, 0, 0; -t71, 0, 0, 0, 0, 0, 0, 0, 0, 0; t70, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiRD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:52
	% EndTime: 2020-05-07 04:44:52
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (178->16), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t85 = -qJ(2) - qJ(7) + pkin(15) - qJ(8);
	t84 = cos(t85);
	t86 = -qJD(2) - qJD(7) - qJD(8);
	t93 = t86 * t84;
	t87 = sin(qJ(1));
	t92 = t86 * t87;
	t88 = cos(qJ(1));
	t91 = t86 * t88;
	t90 = qJD(1) * t87;
	t89 = qJD(1) * t88;
	t83 = sin(t85);
	t82 = t86 * t83;
	t81 = -t83 * t92 + t84 * t89;
	t80 = t83 * t89 + t84 * t92;
	t79 = t83 * t91 + t84 * t90;
	t78 = t83 * t90 - t84 * t91;
	t1 = [t81, t78, 0, 0, 0, 0, t78, t78, 0, 0; t79, -t80, 0, 0, 0, 0, -t80, -t80, 0, 0; 0, t82, 0, 0, 0, 0, t82, t82, 0, 0; t80, -t79, 0, 0, 0, 0, -t79, -t79, 0, 0; t78, t81, 0, 0, 0, 0, t81, t81, 0, 0; 0, -t93, 0, 0, 0, 0, -t93, -t93, 0, 0; -t90, 0, 0, 0, 0, 0, 0, 0, 0, 0; t89, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end