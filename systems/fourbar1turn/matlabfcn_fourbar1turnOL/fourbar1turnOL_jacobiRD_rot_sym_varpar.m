% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% fourbar1turnOL
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = fourbar1turnOL_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1turnOL_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_jacobiRD_rot_sym_varpar: pkin has to be [5x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:41:19
	% EndTime: 2020-04-12 19:41:19
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:41:19
	% EndTime: 2020-04-12 19:41:19
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
	% StartTime: 2020-04-12 19:41:19
	% EndTime: 2020-04-12 19:41:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t31 = sin(qJ(1));
	t38 = qJD(1) * t31;
	t33 = cos(qJ(1));
	t37 = qJD(1) * t33;
	t30 = sin(qJ(2));
	t36 = qJD(2) * t30;
	t32 = cos(qJ(2));
	t35 = qJD(2) * t32;
	t34 = qJD(2) * t33;
	t29 = t31 * t36 - t32 * t37;
	t28 = t30 * t37 + t31 * t35;
	t27 = t30 * t34 + t32 * t38;
	t26 = t30 * t38 - t32 * t34;
	t1 = [t29, t26, 0, 0, 0; -t27, -t28, 0, 0, 0; 0, -t36, 0, 0, 0; t28, t27, 0, 0, 0; t26, t29, 0, 0, 0; 0, -t35, 0, 0, 0; -t38, 0, 0, 0, 0; t37, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:41:19
	% EndTime: 2020-04-12 19:41:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (57->10), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t75 = qJD(2) + qJD(3);
	t77 = sin(qJ(1));
	t82 = t75 * t77;
	t78 = cos(qJ(1));
	t81 = t75 * t78;
	t80 = qJD(1) * t77;
	t79 = qJD(1) * t78;
	t76 = qJ(2) + qJ(3);
	t74 = cos(t76);
	t73 = sin(t76);
	t72 = t75 * t74;
	t71 = t75 * t73;
	t70 = -t73 * t82 + t74 * t79;
	t69 = t73 * t79 + t74 * t82;
	t68 = t73 * t81 + t74 * t80;
	t67 = -t73 * t80 + t74 * t81;
	t1 = [t70, t67, t67, 0, 0; t68, t69, t69, 0, 0; 0, t71, t71, 0, 0; -t69, -t68, -t68, 0, 0; t67, t70, t70, 0, 0; 0, t72, t72, 0, 0; -t80, 0, 0, 0, 0; t79, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:41:19
	% EndTime: 2020-04-12 19:41:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t32 = sin(qJ(1));
	t39 = qJD(1) * t32;
	t34 = cos(qJ(1));
	t38 = qJD(1) * t34;
	t31 = sin(qJ(4));
	t37 = qJD(4) * t31;
	t33 = cos(qJ(4));
	t36 = qJD(4) * t33;
	t35 = qJD(4) * t34;
	t30 = t32 * t37 - t33 * t38;
	t29 = t31 * t38 + t32 * t36;
	t28 = t31 * t35 + t33 * t39;
	t27 = t31 * t39 - t33 * t35;
	t1 = [t30, 0, 0, t27, 0; -t28, 0, 0, -t29, 0; 0, 0, 0, -t37, 0; t29, 0, 0, t28, 0; t27, 0, 0, t30, 0; 0, 0, 0, -t36, 0; -t39, 0, 0, 0, 0; t38, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end