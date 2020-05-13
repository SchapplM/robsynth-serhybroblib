% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh1m2OL
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% 
% Output:
% JRD_rot [9x13]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = palh1m2OL_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),uint8(0),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_jacobiRD_rot_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m2OL_jacobiRD_rot_sym_varpar: qJD has to be [13x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m2OL_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_jacobiRD_rot_sym_varpar: pkin has to be [20x1] (double)');
JRD_rot=NaN(9,13);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:41
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:41
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:41
	% EndTime: 2020-05-02 23:30:42
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->8), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t32 = sin(qJ(1));
	t39 = qJD(1) * t32;
	t34 = cos(qJ(1));
	t38 = qJD(1) * t34;
	t31 = sin(qJ(2));
	t37 = qJD(2) * t31;
	t33 = cos(qJ(2));
	t36 = qJD(2) * t33;
	t35 = qJD(2) * t34;
	t30 = -t32 * t37 + t33 * t38;
	t29 = t31 * t38 + t32 * t36;
	t28 = t31 * t35 + t33 * t39;
	t27 = t31 * t39 - t33 * t35;
	t1 = [t29, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t27, -t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t30, -t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t28, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:42
	% EndTime: 2020-05-02 23:30:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (61->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t68 = qJ(2) + qJ(3);
	t65 = sin(t68);
	t67 = qJD(2) + qJD(3);
	t76 = t67 * t65;
	t66 = cos(t68);
	t75 = t67 * t66;
	t69 = sin(qJ(1));
	t74 = t67 * t69;
	t70 = cos(qJ(1));
	t73 = t67 * t70;
	t72 = qJD(1) * t69;
	t71 = qJD(1) * t70;
	t64 = t65 * t74 - t66 * t71;
	t63 = t65 * t71 + t66 * t74;
	t62 = t65 * t73 + t66 * t72;
	t61 = t65 * t72 - t66 * t73;
	t1 = [t64, t61, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t62, -t63, -t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t76, -t76, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t63, t62, t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t61, t64, t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t75, -t75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t72, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t71, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:42
	% EndTime: 2020-05-02 23:30:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (143->17), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t87 = qJ(2) + qJ(3) + qJ(4);
	t84 = sin(t87);
	t86 = qJD(2) + qJD(3) + qJD(4);
	t95 = t86 * t84;
	t85 = cos(t87);
	t94 = t86 * t85;
	t88 = sin(qJ(1));
	t93 = t86 * t88;
	t89 = cos(qJ(1));
	t92 = t86 * t89;
	t91 = qJD(1) * t88;
	t90 = qJD(1) * t89;
	t83 = t84 * t93 - t85 * t90;
	t82 = t84 * t90 + t85 * t93;
	t81 = t84 * t92 + t85 * t91;
	t80 = t84 * t91 - t85 * t92;
	t1 = [t83, t80, t80, t80, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t81, -t82, -t82, -t82, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t95, -t95, -t95, 0, 0, 0, 0, 0, 0, 0, 0, 0; t82, t81, t81, t81, 0, 0, 0, 0, 0, 0, 0, 0, 0; t80, t83, t83, t83, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t94, -t94, -t94, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t91, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t90, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:43
	% EndTime: 2020-05-02 23:30:44
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (340->29), mult. (279->57), div. (0->0), fcn. (279->6), ass. (0->41)
	t348 = qJD(2) + qJD(3) + qJD(4);
	t350 = sin(qJ(5));
	t372 = t348 * t350;
	t351 = sin(qJ(1));
	t371 = t348 * t351;
	t352 = cos(qJ(5));
	t370 = t348 * t352;
	t353 = cos(qJ(1));
	t369 = t348 * t353;
	t368 = t352 * t353;
	t367 = qJD(1) * t351;
	t366 = qJD(1) * t353;
	t365 = qJD(5) * t350;
	t364 = qJD(5) * t352;
	t363 = qJD(5) * t353;
	t349 = qJ(2) + qJ(3) + qJ(4);
	t346 = sin(t349);
	t362 = t346 * t370;
	t361 = t346 * t371;
	t347 = cos(t349);
	t360 = t347 * t371;
	t359 = t346 * t369;
	t358 = t347 * t369;
	t357 = qJD(5) * t347 - qJD(1);
	t356 = qJD(1) * t347 - qJD(5);
	t355 = t357 * t350;
	t354 = t356 * t351 + t359;
	t345 = t348 * t347;
	t344 = t347 * t366 - t361;
	t343 = -t347 * t367 - t359;
	t342 = -t347 * t365 - t362;
	t341 = t346 * t372 - t347 * t364;
	t340 = -t352 * t360 + (t351 * t365 - t352 * t366) * t346;
	t339 = t350 * t360 + (t350 * t366 + t351 * t364) * t346;
	t338 = -t352 * t358 + (t350 * t363 + t352 * t367) * t346;
	t337 = t350 * t358 + (-t350 * t367 + t352 * t363) * t346;
	t336 = -t356 * t368 + (t355 + t362) * t351;
	t335 = t357 * t352 * t351 + (t356 * t353 - t361) * t350;
	t334 = t354 * t352 + t353 * t355;
	t333 = t354 * t350 - t357 * t368;
	t1 = [t336, t338, t338, t338, t333, 0, 0, 0, 0, 0, 0, 0, 0; -t334, t340, t340, t340, -t335, 0, 0, 0, 0, 0, 0, 0, 0; 0, t342, t342, t342, -t346 * t364 - t347 * t372, 0, 0, 0, 0, 0, 0, 0, 0; t335, t337, t337, t337, t334, 0, 0, 0, 0, 0, 0, 0, 0; t333, t339, t339, t339, t336, 0, 0, 0, 0, 0, 0, 0, 0; 0, t341, t341, t341, t346 * t365 - t347 * t370, 0, 0, 0, 0, 0, 0, 0, 0; -t346 * t366 - t360, t343, t343, t343, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t346 * t367 + t358, t344, t344, t344, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t345, t345, t345, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:42
	% EndTime: 2020-05-02 23:30:42
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
	t1 = [t30, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0; -t28, 0, 0, 0, 0, -t29, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t37, 0, 0, 0, 0, 0, 0, 0; t29, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0; t27, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t36, 0, 0, 0, 0, 0, 0, 0; -t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiRD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:42
	% EndTime: 2020-05-02 23:30:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (59->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t77 = qJ(2) + qJ(7);
	t75 = cos(t77);
	t76 = qJD(2) + qJD(7);
	t84 = t76 * t75;
	t78 = sin(qJ(1));
	t83 = t76 * t78;
	t79 = cos(qJ(1));
	t82 = t76 * t79;
	t81 = qJD(1) * t78;
	t80 = qJD(1) * t79;
	t74 = sin(t77);
	t73 = t76 * t74;
	t70 = -t74 * t83 + t75 * t80;
	t69 = t74 * t80 + t75 * t83;
	t68 = t74 * t82 + t75 * t81;
	t67 = t74 * t81 - t75 * t82;
	t1 = [t69, t68, 0, 0, 0, 0, t68, 0, 0, 0, 0, 0, 0; t67, -t70, 0, 0, 0, 0, -t70, 0, 0, 0, 0, 0, 0; 0, -t84, 0, 0, 0, 0, -t84, 0, 0, 0, 0, 0, 0; t70, -t67, 0, 0, 0, 0, -t67, 0, 0, 0, 0, 0, 0; t68, t69, 0, 0, 0, 0, t69, 0, 0, 0, 0, 0, 0; 0, t73, 0, 0, 0, 0, t73, 0, 0, 0, 0, 0, 0; -t81, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiRD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:42
	% EndTime: 2020-05-02 23:30:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (59->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t49 = qJ(2) + qJ(8);
	t47 = cos(t49);
	t48 = qJD(2) + qJD(8);
	t56 = t48 * t47;
	t50 = sin(qJ(1));
	t55 = t48 * t50;
	t51 = cos(qJ(1));
	t54 = t48 * t51;
	t53 = qJD(1) * t50;
	t52 = qJD(1) * t51;
	t46 = sin(t49);
	t45 = t48 * t46;
	t42 = -t46 * t55 + t47 * t52;
	t41 = t46 * t52 + t47 * t55;
	t40 = t46 * t54 + t47 * t53;
	t39 = t46 * t53 - t47 * t54;
	t1 = [t41, t40, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0; t39, -t42, 0, 0, 0, 0, 0, -t42, 0, 0, 0, 0, 0; 0, -t56, 0, 0, 0, 0, 0, -t56, 0, 0, 0, 0, 0; t42, -t39, 0, 0, 0, 0, 0, -t39, 0, 0, 0, 0, 0; t40, t41, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0; 0, t45, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0; -t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiRD_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:42
	% EndTime: 2020-05-02 23:30:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (140->12), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t79 = qJ(2) + qJ(8) + qJ(9);
	t76 = sin(t79);
	t78 = qJD(2) + qJD(8) + qJD(9);
	t86 = t78 * t76;
	t80 = sin(qJ(1));
	t85 = t78 * t80;
	t81 = cos(qJ(1));
	t84 = t78 * t81;
	t83 = qJD(1) * t80;
	t82 = qJD(1) * t81;
	t77 = cos(t79);
	t75 = t78 * t77;
	t74 = -t76 * t85 + t77 * t82;
	t73 = -t76 * t82 - t77 * t85;
	t72 = -t76 * t84 - t77 * t83;
	t71 = t76 * t83 - t77 * t84;
	t1 = [t73, t72, 0, 0, 0, 0, 0, t72, t72, 0, 0, 0, 0; -t71, t74, 0, 0, 0, 0, 0, t74, t74, 0, 0, 0, 0; 0, t75, 0, 0, 0, 0, 0, t75, t75, 0, 0, 0, 0; -t74, t71, 0, 0, 0, 0, 0, t71, t71, 0, 0, 0, 0; t72, t73, 0, 0, 0, 0, 0, t73, t73, 0, 0, 0, 0; 0, -t86, 0, 0, 0, 0, 0, -t86, -t86, 0, 0, 0, 0; -t83, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t82, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiRD_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:42
	% EndTime: 2020-05-02 23:30:42
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (181->17), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t93 = -qJ(2) - qJ(7) + pkin(19) - qJ(10);
	t91 = sin(t93);
	t94 = -qJD(2) - qJD(7) - qJD(10);
	t102 = t94 * t91;
	t92 = cos(t93);
	t101 = t94 * t92;
	t95 = sin(qJ(1));
	t100 = t94 * t95;
	t96 = cos(qJ(1));
	t99 = t94 * t96;
	t98 = qJD(1) * t95;
	t97 = qJD(1) * t96;
	t89 = -t91 * t100 + t92 * t97;
	t88 = t92 * t100 + t91 * t97;
	t87 = -t91 * t99 - t92 * t98;
	t86 = t91 * t98 - t92 * t99;
	t1 = [t88, t87, 0, 0, 0, 0, t87, 0, 0, t87, 0, 0, 0; t86, t89, 0, 0, 0, 0, t89, 0, 0, t89, 0, 0, 0; 0, -t101, 0, 0, 0, 0, -t101, 0, 0, -t101, 0, 0, 0; -t89, -t86, 0, 0, 0, 0, -t86, 0, 0, -t86, 0, 0, 0; t87, t88, 0, 0, 0, 0, t88, 0, 0, t88, 0, 0, 0; 0, -t102, 0, 0, 0, 0, -t102, 0, 0, -t102, 0, 0, 0; -t98, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end