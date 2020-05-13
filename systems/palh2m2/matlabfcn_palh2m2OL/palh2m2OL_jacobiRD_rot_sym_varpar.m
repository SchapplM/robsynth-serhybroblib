% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = palh2m2OL_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m2OL_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_jacobiRD_rot_sym_varpar: pkin has to be [5x1] (double)');
JRD_rot=NaN(9,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
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
	t1 = [t30, t27, 0, 0, 0, 0; -t28, -t29, 0, 0, 0, 0; 0, -t37, 0, 0, 0, 0; t29, t28, 0, 0, 0, 0; t27, t30, 0, 0, 0, 0; 0, -t36, 0, 0, 0, 0; -t39, 0, 0, 0, 0, 0; t38, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (61->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t67 = qJ(2) + qJ(3);
	t64 = sin(t67);
	t66 = qJD(2) + qJD(3);
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
	t1 = [t63, t60, t60, 0, 0, 0; -t61, -t62, -t62, 0, 0, 0; 0, -t75, -t75, 0, 0, 0; t62, t61, t61, 0, 0, 0; t60, t63, t63, 0, 0, 0; 0, -t74, -t74, 0, 0, 0; -t71, 0, 0, 0, 0, 0; t70, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:58
	% EndTime: 2020-05-03 06:34:58
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (143->17), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t86 = qJ(2) + qJ(3) + qJ(4);
	t83 = sin(t86);
	t85 = qJD(2) + qJD(3) + qJD(4);
	t94 = t85 * t83;
	t84 = cos(t86);
	t93 = t85 * t84;
	t87 = sin(qJ(1));
	t92 = t85 * t87;
	t88 = cos(qJ(1));
	t91 = t85 * t88;
	t90 = qJD(1) * t87;
	t89 = qJD(1) * t88;
	t82 = t83 * t92 - t84 * t89;
	t81 = t83 * t89 + t84 * t92;
	t80 = t83 * t91 + t84 * t90;
	t79 = t83 * t90 - t84 * t91;
	t1 = [t82, t79, t79, t79, 0, 0; -t80, -t81, -t81, -t81, 0, 0; 0, -t94, -t94, -t94, 0, 0; t81, t80, t80, t80, 0, 0; t79, t82, t82, t82, 0, 0; 0, -t93, -t93, -t93, 0, 0; -t90, 0, 0, 0, 0, 0; t89, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:58
	% EndTime: 2020-05-03 06:34:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (257->20), mult. (90->14), div. (0->0), fcn. (90->4), ass. (0->17)
	t101 = qJD(2) + qJD(3) + qJD(4) + qJD(5);
	t102 = qJ(2) + qJ(3) + qJ(4) + qJ(5);
	t99 = sin(t102);
	t110 = t101 * t99;
	t100 = cos(t102);
	t109 = t101 * t100;
	t103 = sin(qJ(1));
	t108 = t101 * t103;
	t104 = cos(qJ(1));
	t107 = t101 * t104;
	t106 = qJD(1) * t103;
	t105 = qJD(1) * t104;
	t98 = -t100 * t105 + t99 * t108;
	t97 = t100 * t108 + t99 * t105;
	t96 = t100 * t106 + t99 * t107;
	t95 = -t100 * t107 + t99 * t106;
	t1 = [t98, t95, t95, t95, t95, 0; -t96, -t97, -t97, -t97, -t97, 0; 0, -t110, -t110, -t110, -t110, 0; t97, t96, t96, t96, t96, 0; t95, t98, t98, t98, t98, 0; 0, -t109, -t109, -t109, -t109, 0; -t106, 0, 0, 0, 0, 0; t105, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:59
	% EndTime: 2020-05-03 06:35:00
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (580->33), mult. (332->57), div. (0->0), fcn. (332->6), ass. (0->41)
	t337 = qJ(2) + qJ(3) + qJ(4) + qJ(5);
	t335 = cos(t337);
	t336 = qJD(2) + qJD(3) + qJD(4) + qJD(5);
	t361 = t336 * t335;
	t338 = sin(qJ(6));
	t360 = t336 * t338;
	t339 = sin(qJ(1));
	t359 = t336 * t339;
	t340 = cos(qJ(6));
	t358 = t336 * t340;
	t341 = cos(qJ(1));
	t357 = t336 * t341;
	t356 = t340 * t341;
	t355 = qJD(1) * t339;
	t354 = qJD(1) * t341;
	t353 = qJD(6) * t338;
	t352 = qJD(6) * t340;
	t351 = qJD(6) * t341;
	t334 = sin(t337);
	t350 = t334 * t358;
	t349 = t334 * t359;
	t348 = t335 * t359;
	t347 = t334 * t357;
	t346 = t335 * t357;
	t345 = qJD(6) * t335 + qJD(1);
	t344 = qJD(1) * t335 + qJD(6);
	t343 = t345 * t338;
	t342 = t344 * t339 + t347;
	t333 = -t335 * t354 + t349;
	t332 = t335 * t355 + t347;
	t331 = -t335 * t353 - t350;
	t330 = t334 * t360 - t335 * t352;
	t329 = -t340 * t348 + (t339 * t353 - t340 * t354) * t334;
	t328 = t338 * t348 + (t338 * t354 + t339 * t352) * t334;
	t327 = -t340 * t346 + (t338 * t351 + t340 * t355) * t334;
	t326 = t338 * t346 + (-t338 * t355 + t340 * t351) * t334;
	t325 = -t344 * t356 + (t343 + t350) * t339;
	t324 = t345 * t340 * t339 + (t344 * t341 - t349) * t338;
	t323 = t342 * t340 + t341 * t343;
	t322 = t342 * t338 - t345 * t356;
	t1 = [t325, t327, t327, t327, t327, t322; -t323, t329, t329, t329, t329, -t324; 0, t331, t331, t331, t331, -t334 * t352 - t335 * t360; t324, t326, t326, t326, t326, t323; t322, t328, t328, t328, t328, t325; 0, t330, t330, t330, t330, t334 * t353 - t335 * t358; t334 * t354 + t348, t332, t332, t332, t332, 0; t334 * t355 - t346, t333, t333, t333, t333, 0; 0, -t361, -t361, -t361, -t361, 0;];
	JRD_rot = t1;
end