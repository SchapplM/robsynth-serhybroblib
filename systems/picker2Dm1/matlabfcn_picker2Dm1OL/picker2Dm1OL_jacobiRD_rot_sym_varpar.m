% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% picker2Dm1OL
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% 
% Output:
% JRD_rot [9x12]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = picker2Dm1OL_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_jacobiRD_rot_sym_varpar: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm1OL_jacobiRD_rot_sym_varpar: qJD has to be [12x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'picker2Dm1OL_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
JRD_rot=NaN(9,12);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-11 05:46:38
	% EndTime: 2020-05-11 05:46:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-11 05:46:38
	% EndTime: 2020-05-11 05:46:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t30 = qJD(1) * sin(qJ(1));
	t28 = qJD(1) * cos(qJ(1));
	t1 = [t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-11 05:46:38
	% EndTime: 2020-05-11 05:46:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->4), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t54 = qJ(1) + qJ(2);
	t53 = qJD(1) + qJD(2);
	t51 = t53 * cos(t54);
	t50 = t53 * sin(t54);
	t1 = [t51, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t50, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t50, -t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t51, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-11 05:46:38
	% EndTime: 2020-05-11 05:46:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (57->11), mult. (12->2), div. (0->0), fcn. (12->2), ass. (0->5)
	t59 = qJD(1) + qJD(2) + qJD(3);
	t60 = qJ(1) + qJ(2) + qJ(3);
	t61 = t59 * cos(t60);
	t56 = t59 * sin(t60);
	t1 = [-t61, -t61, -t61, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t56, -t56, -t56, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t56, t56, t56, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t61, -t61, -t61, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-11 05:46:38
	% EndTime: 2020-05-11 05:46:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (51->5), mult. (12->2), div. (0->0), fcn. (12->2), ass. (0->5)
	t66 = qJ(1) + qJ(2) + qJ(4);
	t65 = qJD(1) + qJD(2) + qJD(4);
	t63 = t65 * cos(t66);
	t62 = t65 * sin(t66);
	t1 = [t63, t63, 0, t63, 0, 0, 0, 0, 0, 0, 0, 0; t62, t62, 0, t62, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t62, -t62, 0, -t62, 0, 0, 0, 0, 0, 0, 0, 0; t63, t63, 0, t63, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-11 05:46:38
	% EndTime: 2020-05-11 05:46:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = pkin(8) + qJ(5);
	t37 = qJD(5) * sin(t35);
	t36 = qJD(5) * cos(t35);
	t1 = [0, 0, 0, 0, -t36, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t37, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t36, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-11 05:46:38
	% EndTime: 2020-05-11 05:46:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (57->11), mult. (12->2), div. (0->0), fcn. (12->2), ass. (0->5)
	t47 = qJD(1) + qJD(2) + qJD(6);
	t48 = qJ(1) + qJ(2) + qJ(6);
	t49 = t47 * cos(t48);
	t44 = t47 * sin(t48);
	t1 = [-t49, -t49, 0, 0, 0, -t49, 0, 0, 0, 0, 0, 0; -t44, -t44, 0, 0, 0, -t44, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t44, t44, 0, 0, 0, t44, 0, 0, 0, 0, 0, 0; -t49, -t49, 0, 0, 0, -t49, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiRD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-11 05:46:38
	% EndTime: 2020-05-11 05:46:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(7) * sin(qJ(7));
	t30 = qJD(7) * cos(qJ(7));
	t1 = [0, 0, 0, 0, 0, 0, -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiRD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-11 05:46:38
	% EndTime: 2020-05-11 05:46:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t47 = qJD(1) + qJD(8);
	t48 = qJ(1) + qJ(8);
	t49 = t47 * cos(t48);
	t44 = t47 * sin(t48);
	t1 = [-t49, 0, 0, 0, 0, 0, 0, -t49, 0, 0, 0, 0; -t44, 0, 0, 0, 0, 0, 0, -t44, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t44, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0; -t49, 0, 0, 0, 0, 0, 0, -t49, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiRD_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-11 05:46:38
	% EndTime: 2020-05-11 05:46:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (100->6), mult. (16->2), div. (0->0), fcn. (16->2), ass. (0->5)
	t78 = qJ(1) + qJ(2) + qJ(3) + qJ(9);
	t77 = qJD(1) + qJD(2) + qJD(3) + qJD(9);
	t75 = t77 * cos(t78);
	t74 = t77 * sin(t78);
	t1 = [t75, t75, t75, 0, 0, 0, 0, 0, t75, 0, 0, 0; t74, t74, t74, 0, 0, 0, 0, 0, t74, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t74, -t74, -t74, 0, 0, 0, 0, 0, -t74, 0, 0, 0; t75, t75, t75, 0, 0, 0, 0, 0, t75, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiRD_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-11 05:46:38
	% EndTime: 2020-05-11 05:46:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (108->14), mult. (16->2), div. (0->0), fcn. (16->2), ass. (0->5)
	t71 = qJD(1) + qJD(2) + qJD(4) + qJD(10);
	t72 = qJ(1) + qJ(2) + qJ(4) + qJ(10);
	t73 = t71 * cos(t72);
	t68 = t71 * sin(t72);
	t1 = [-t73, -t73, 0, -t73, 0, 0, 0, 0, 0, -t73, 0, 0; -t68, -t68, 0, -t68, 0, 0, 0, 0, 0, -t68, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t68, t68, 0, t68, 0, 0, 0, 0, 0, t68, 0, 0; -t73, -t73, 0, -t73, 0, 0, 0, 0, 0, -t73, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end