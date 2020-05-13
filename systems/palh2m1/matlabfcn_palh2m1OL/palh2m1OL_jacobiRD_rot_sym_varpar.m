% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = palh2m1OL_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m1OL_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_jacobiRD_rot_sym_varpar: pkin has to be [6x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->7), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
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
	t1 = [t30, t27, 0, 0, 0; -t28, -t29, 0, 0, 0; 0, t37, 0, 0, 0; t29, t28, 0, 0, 0; t27, t30, 0, 0, 0; 0, t36, 0, 0, 0; t39, 0, 0, 0, 0; -t38, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (57->10), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t72 = qJD(2) + qJD(3);
	t74 = sin(qJ(1));
	t79 = t72 * t74;
	t75 = cos(qJ(1));
	t78 = t72 * t75;
	t77 = qJD(1) * t74;
	t76 = qJD(1) * t75;
	t73 = qJ(2) + qJ(3);
	t71 = cos(t73);
	t70 = sin(t73);
	t69 = t72 * t71;
	t68 = t72 * t70;
	t67 = t70 * t79 - t71 * t76;
	t66 = t70 * t76 + t71 * t79;
	t65 = t70 * t78 + t71 * t77;
	t64 = t70 * t77 - t71 * t78;
	t1 = [t67, t64, t64, 0, 0; -t65, -t66, -t66, 0, 0; 0, t68, t68, 0, 0; t66, t65, t65, 0, 0; t64, t67, t67, 0, 0; 0, t69, t69, 0, 0; t77, 0, 0, 0, 0; -t76, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (137->11), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t91 = qJD(2) + qJD(3) + qJD(4);
	t93 = sin(qJ(1));
	t98 = t91 * t93;
	t94 = cos(qJ(1));
	t97 = t91 * t94;
	t96 = qJD(1) * t93;
	t95 = qJD(1) * t94;
	t92 = qJ(2) + qJ(3) + qJ(4);
	t90 = cos(t92);
	t89 = sin(t92);
	t88 = t91 * t90;
	t87 = t91 * t89;
	t86 = t89 * t98 - t90 * t95;
	t85 = t89 * t95 + t90 * t98;
	t84 = t89 * t97 + t90 * t96;
	t83 = t89 * t96 - t90 * t97;
	t1 = [t86, t83, t83, t83, 0; -t84, -t85, -t85, -t85, 0; 0, t87, t87, t87, 0; t85, t84, t84, t84, 0; t83, t86, t86, t86, 0; 0, t88, t88, t88, 0; t96, 0, 0, 0, 0; -t95, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:43
	% EndTime: 2020-05-03 00:53:43
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (343->32), mult. (279->57), div. (0->0), fcn. (279->6), ass. (0->41)
	t343 = qJ(2) + qJ(3) + qJ(4);
	t341 = cos(t343);
	t342 = qJD(2) + qJD(3) + qJD(4);
	t367 = t342 * t341;
	t344 = sin(qJ(5));
	t366 = t342 * t344;
	t345 = sin(qJ(1));
	t365 = t342 * t345;
	t346 = cos(qJ(5));
	t364 = t342 * t346;
	t347 = cos(qJ(1));
	t363 = t342 * t347;
	t362 = t346 * t347;
	t361 = qJD(1) * t345;
	t360 = qJD(1) * t347;
	t359 = qJD(5) * t344;
	t358 = qJD(5) * t346;
	t357 = qJD(5) * t347;
	t340 = sin(t343);
	t356 = t340 * t364;
	t355 = t340 * t365;
	t354 = t341 * t365;
	t353 = t340 * t363;
	t352 = t341 * t363;
	t351 = qJD(5) * t341 + qJD(1);
	t350 = qJD(1) * t341 + qJD(5);
	t349 = t351 * t344;
	t348 = t350 * t345 + t353;
	t339 = t341 * t360 - t355;
	t338 = -t341 * t361 - t353;
	t337 = t341 * t359 + t356;
	t336 = -t340 * t366 + t341 * t358;
	t335 = -t346 * t354 + (t345 * t359 - t346 * t360) * t340;
	t334 = t344 * t354 + (t344 * t360 + t345 * t358) * t340;
	t333 = -t346 * t352 + (t344 * t357 + t346 * t361) * t340;
	t332 = t344 * t352 + (-t344 * t361 + t346 * t357) * t340;
	t331 = -t350 * t362 + (t349 + t356) * t345;
	t330 = t351 * t346 * t345 + (t350 * t347 - t355) * t344;
	t329 = t348 * t346 + t347 * t349;
	t328 = t348 * t344 - t351 * t362;
	t1 = [t331, t333, t333, t333, t328; -t329, t335, t335, t335, -t330; 0, t337, t337, t337, t340 * t358 + t341 * t366; t330, t332, t332, t332, t329; t328, t334, t334, t334, t331; 0, t336, t336, t336, -t340 * t359 + t341 * t364; -t340 * t360 - t354, t338, t338, t338, 0; -t340 * t361 + t352, t339, t339, t339, 0; 0, -t367, -t367, -t367, 0;];
	JRD_rot = t1;
end